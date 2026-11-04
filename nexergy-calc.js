/* Nexergy 계산 파일 — nexergy_engine.py (v0.1) 의 열저항망 모델을 그대로 옮긴 것.
   외기 → 외관(재킷) → [복사 ‖ 잔류가스 전도 ‖ 지지대 전도] → 내관 벽 → 유체
   Q_total = Σ 경로 W + 부품 W,  BOG[kg/년] = Q × 운전시간 / 증발잠열,  손실 = BOG × 단가 × (1−회수율)

   화면(app.js)이 넘겨주는 입력 p 와 돌려주는 결과의 이름은 '계산 연결 규칙'을 따른다.
   파이썬 엔진과 다른 점은 없고, 화면에 없는 입력(외관 치수·방사율·열화계수 등)은
   아래 ADV 기본값을 쓴다. 팀에서 확정되면 ADV 만 바꾸면 된다. */
window.NexergyCalc = (function () {
  'use strict';

  var SIGMA = 5.670374419e-8;   // W/m²K⁴
  var R_GAS = 8.314462618;      // J/mol·K

  /* ---- 화면에 없는 고급 입력 (nexergy_engine.py 의 Insulation/Operating/Segment 기본값) ---- */
  var ADV = {
    h_external: 10.0,        // 외관 바깥 자연대류+복사 등가 열전달계수 W/m²K
    emissivity_mli: 0.03,    // 알루미늄 증착 필름 방사율
    emissivity_pipe: 0.10,   // 스테인리스 표면 방사율 (진공만)
    mli_degradation: 3.0,    // MLI 열화계수 D (실측 문헌 2~5)
    mli_layer_pitch_mm: 0.5, // 층 간격 (연속체 가스전도 두께용)
    accommodation: 0.9,      // 열적응계수 α
    gas_molar_mass: 0.029,   // 잔류가스 몰질량 kg/mol (공기)
    gas_gamma: 1.4,          // 비열비
    gas_gauge_temp_k: 300.0, // 압력계 기준 온도
    jacket_gap_mm: 60.0,     // 화면에 외관 외경이 없을 때: 외관 외경 = 내관 외경 + 60 mm
    jacket_wt_mm: 3.0,       // 화면에 외관 두께가 없을 때
    inner_wt_mm: 2.0,        // 화면에 내관 두께가 없을 때
    support_ref_dt: 278.0    // NIST 식이 없는 지지대 재질: 화면 kint 값을 이 온도차 기준으로 비례
  };

  /* ---- 1. 유체 물성: NIST 몇 점 보간 ---- */
  var TSAT = {   // (절대압 bar, 포화온도 K)
    LH2: [[1.013, 20.3], [2.0, 22.8], [3.0, 24.6], [5.0, 27.1], [8.0, 29.9], [10.0, 31.4]],
    LN2: [[1.013, 77.4], [2.0, 84.0], [3.0, 88.0], [5.0, 94.0], [8.0, 100.4], [10.0, 103.7]],
    LNG: [[1.013, 111.7], [2.0, 120.6], [3.0, 126.5], [5.0, 134.7], [8.0, 143.6], [10.0, 149.1]]
  };
  var HFG = {    // (포화온도 K, 증발잠열 kJ/kg)
    LH2: [[20.3, 446.0], [23.0, 430.0], [25.0, 412.0], [27.0, 390.0], [30.0, 335.0], [32.0, 260.0]],
    LN2: [[77.4, 199.0], [85.0, 187.0], [95.0, 170.0], [105.0, 150.0]],
    LNG: [[111.7, 510.0], [120.0, 495.0], [135.0, 465.0], [150.0, 425.0]]
  };

  function interp(x, pts, logx) {
    var xs = pts.map(function (q) { return logx ? Math.log(q[0]) : q[0]; });
    var ys = pts.map(function (q) { return q[1]; });
    if (logx) x = Math.log(x);
    if (x <= xs[0]) return ys[0];
    if (x >= xs[xs.length - 1]) return ys[ys.length - 1];
    for (var i = 0; i < xs.length - 1; i++) {
      if (xs[i] <= x && x <= xs[i + 1]) {
        var t = (x - xs[i]) / (xs[i + 1] - xs[i]);
        return ys[i] + t * (ys[i + 1] - ys[i]);
      }
    }
    return ys[ys.length - 1];
  }
  function saturationTemperature(fluid, pBarg) { return interp(pBarg + 1.013, TSAT[fluid], true); }   // K
  function latentHeat(fluid, tSat) { return interp(tSat, HFG[fluid], false) * 1000.0; }              // J/kg

  /* ---- 2. 고체 열전도도: NIST 극저온 물성 적합식 log10(k) = Σ a_i (log10 T)^i ---- */
  var NIST_FIT = {
    G10:   [-4.1236, 13.788, -26.068, 26.272, -14.663, 4.4954, -0.6905, 0.0397, 0.0],
    PTFE:  [2.7380, -30.677, 89.430, -136.99, 124.69, -69.556, 23.320, -4.3135, 0.33829],
    SS304: [-1.4087, 1.3982, 0.2543, -0.6260, 0.2334, 0.4256, -0.4658, 0.1650, -0.0199]
  };
  var MAT_ALIAS = { G10: 'G10', 'G-10': 'G10', PTFE: 'PTFE', SUS304: 'SS304', SS304: 'SS304', SUS316: 'SS304', SUS316L: 'SS304', SS316: 'SS304' };

  function conductivity(mat, T) {
    var c = NIST_FIT[mat];
    T = Math.min(Math.max(T, 4.0), 300.0);
    var lt = Math.log10(T), s = 0;
    for (var i = 0; i < c.length; i++) s += c[i] * Math.pow(lt, i);
    return Math.pow(10, s);
  }
  function conductivityIntegral(mat, tc, th, n) {   // ∫k dT [W/m], 사다리꼴
    n = n || 200;
    var h = (th - tc) / n, s = 0.5 * (conductivity(mat, tc) + conductivity(mat, th));
    for (var i = 1; i < n; i++) s += conductivity(mat, tc + i * h);
    return s * h;
  }
  function gasConductivityAir(T) { return 0.0241 * Math.pow(T / 273.15, 0.9); }

  /* ---- 3. 구간 1개 열유입 (nexergy_engine._segment_heat) ---- */
  function segmentHeat(s, p, tCold, tAmb, supKintFn) {
    var L = s.L;
    var innerOd = s.D * 1000.0;                                           // mm
    var outerOd = s.ood > 0 ? s.ood : innerOd + ADV.jacket_gap_mm;
    var outerWt = s.ot > 0 ? s.ot : ADV.jacket_wt_mm;
    var innerWt = s.it > 0 ? s.it : ADV.inner_wt_mm;
    var rI = innerOd / 2000.0;                                            // 내관 외반경 m
    var rJin = (outerOd / 2.0 - outerWt) / 1000.0;
    var rJout = outerOd / 2000.0;
    if (!(rJin > rI)) rJin = rI + ADV.jacket_gap_mm / 2000.0;            // 치수 오류 방어
    var Ai = 2 * Math.PI * rI * L, AJin = 2 * Math.PI * rJin * L, Aout = 2 * Math.PI * rJout * L;

    var nSup = p.spacing > 0 ? Math.max(Math.floor(L / p.spacing + 1e-9), 1) : 0;
    var tHot = tAmb, dT, qRad = 0, qGas = 0, qSup = 0, qFoam = 0, qIns = 0, nGaps = 1, Rext;

    for (var it = 0; it < 4; it++) {
      dT = tHot - tCold;
      qRad = qGas = qSup = qFoam = 0;

      if (p.method === 'foam') {
        var rF = rI + p.foamT;
        var Rfoam = Math.log(rF / rI) / (2 * Math.PI * p.foamK * L);
        Rext = 1.0 / (ADV.h_external * 2 * Math.PI * rF * L);
        qFoam = (tAmb - tCold) / (Rfoam + Rext);
        tHot = tAmb - qFoam * Rext;
        qIns = qFoam;
        nGaps = 1;
      } else {
        // (a) 복사
        if (p.method === 'mli') {
          var N = Math.max(Math.round(p.N), 0), e = ADV.emissivity_mli;
          var epsEff = ADV.mli_degradation * e / ((N + 1) * (2.0 - e));
          qRad = SIGMA * Ai * (Math.pow(tHot, 4) - Math.pow(tCold, 4)) * epsEff;
          nGaps = N + 1;
        } else {   // 진공만: 동심 원통 회색체
          var ep = ADV.emissivity_pipe;
          var denom = 1.0 / ep + (Ai / AJin) * (1.0 / ep - 1.0);
          qRad = SIGMA * Ai * (Math.pow(tHot, 4) - Math.pow(tCold, 4)) / denom;
          nGaps = 1;
        }
        // (b) 잔류가스 전도: 자유분자(Corruccini) ↔ 연속체, 조화평균
        var Lam = (ADV.gas_gamma + 1) / (ADV.gas_gamma - 1) *
                  Math.sqrt(R_GAS / (8 * Math.PI * ADV.gas_molar_mass * ADV.gas_gauge_temp_k));
        var qFm = ADV.accommodation * Lam * p.P * Ai * dT / nGaps;
        var kG = gasConductivityAir(0.5 * (tHot + tCold)), qC;
        if (p.method === 'mli') {
          var tBlanket = Math.max(p.N * ADV.mli_layer_pitch_mm / 1000.0, 1e-4);
          qC = kG * Ai * dT / tBlanket;
        } else {
          qC = 2 * Math.PI * kG * L * dT / Math.log(rJin / rI);
        }
        qGas = qFm <= 0 ? 0 : 1.0 / (1.0 / qFm + 1.0 / qC);
        // (c) 지지대 전도
        if (nSup > 0) qSup = nSup * supKintFn(tCold, tHot) * (p.supArea * 1e-6) / (p.supLen / 1000.0);

        qIns = qRad + qGas + qSup;
        Rext = 1.0 / (ADV.h_external * Aout);
        tHot = tAmb - qIns * Rext;
      }
    }

    // 내관 벽 전도 (참고용, 매우 작음)
    var wallMat = MAT_ALIAS[s.mat] || 'SS304', Rwall = 0;
    if (NIST_FIT[wallMat]) {
      var rInIn = rI - innerWt / 1000.0;
      if (rInIn > 0) Rwall = Math.log(rI / rInIn) / (2 * Math.PI * conductivity(wallMat, tCold + 2.0) * L);
    }

    return {
      name: s.name, tag: s.tag, L: L,
      rad: qRad, gas: qGas, sup: qSup, foam: qFoam, nsup: nSup, Q: qIns,
      perM: qIns / L, tJacket: tHot, dTwall: qIns * Rwall
    };
  }

  /* ---- 4. 전체 계산 (nexergy_engine.calculate) ---- */
  function core(p) {
    var fluid = (p.fluid && p.fluid.short) || 'LH2';
    if (!TSAT[fluid]) fluid = 'LH2';
    var tCold = saturationTemperature(fluid, p.pressure > 0 ? p.pressure : 0);
    var hfg = latentHeat(fluid, tCold);
    var tAmb = p.Ta + 273.15;
    if (!(tAmb > tCold)) tAmb = tCold + 1e-6;

    // 지지대 ∫k dT: NIST 식이 있는 재질은 적분, 없으면 화면의 kint 값을 온도차에 비례
    var supMat = MAT_ALIAS[p.supMat];
    var supKintFn = supMat
      ? function (tc, th) { return conductivityIntegral(supMat, tc, th); }
      : function (tc, th) { return p.supKint * (th - tc) / ADV.support_ref_dt; };

    var segs = p.segs.map(function (s) { return segmentHeat(s, p, tCold, tAmb, supKintFn); });
    var parts = p.parts.map(function (x) { return { type: x.type, qty: x.qty, w: x.w, Q: x.qty * x.w }; });

    function sum(a, k) { return a.reduce(function (t, x) { return t + (+x[k] || 0); }, 0); }
    var L = sum(segs, 'L'), nsup = sum(segs, 'nsup'), partQty = sum(parts, 'qty');
    var Qp = sum(parts, 'Q'), Qr = sum(segs, 'rad'), Qg = sum(segs, 'gas'), Qs = sum(segs, 'sup'), Qf = sum(segs, 'foam');
    var Q = Qp + Qr + Qg + Qs + Qf;

    var secs = p.hours * 3600.0;
    var bogYr = Q * secs / hfg;                 // kg/년 (연간 운전시간 기준)
    var bogDay = Q * 86400.0 / hfg;             // kg/일 (24시간 연속 운전 기준)
    var loss = bogYr * p.price * (1.0 - p.rec / 100.0);

    // 화면 표시용 p: 포화온도·잠열을 운전 압력에서 구한 값으로 덮어씀
    var p2 = Object.assign({}, p, {
      Tf: Math.round((tCold - 273.15) * 10) / 10,
      latent: Math.round(hfg / 1000.0),
      Tsat_K: tCold
    });

    return {
      p: p2, dT: tAmb - tCold,
      segs: segs, parts: parts,
      L: L, nsup: nsup, partQty: partQty,
      Q: Q, perM: L > 0 ? Q / L : 0,
      bogDay: bogDay, bogYr: bogYr, loss: loss,
      Qp: Qp, Qr: Qr, Qg: Qg, Qs: Qs, Qf: Qf
    };
  }

  return {
    name: 'Nexergy 열저항망 엔진 v0.1',
    core: core,
    adv: ADV,
    // 검증·디버그용
    saturationTemperature: saturationTemperature, latentHeat: latentHeat,
    conductivity: conductivity, conductivityIntegral: conductivityIntegral
  };
})();
