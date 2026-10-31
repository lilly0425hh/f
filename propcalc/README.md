# propcalc

고체 로켓 추진제(산화제 + 연료)의 연소실 화학평형, 단열 화염온도, 생성물 조성(응축상 포함), 특성속도 c*, 추력계수 Cf, 비추력 Isp를 계산하는 Python 패키지입니다. 학술 탐구 보고서용 열화학 계산이 목적이며, 추진제 제조나 모터 제작 절차는 범위에 포함하지 않습니다.

## 설치와 실행

Python 3.11 이상이 필요합니다.

```bash
pip install -e ".[cea,rocketcea,dev]"
pytest
python notebooks/stage5_analysis.py
```

- 1순위 백엔드는 NASA 공식 `cea`(Apache-2.0), 교차검증용은 `rocketcea`(GPL-3.0, 소스 빌드에 gfortran 필요)입니다. 이 저장소는 cea 3.3.4, rocketcea 1.2.3으로 검증했습니다.
- `data/ingredients.yaml`은 저장소 체크아웃 기준 경로로 읽습니다. 휠로 설치하면 데이터 파일이 포함되지 않으므로 `load_ingredients(path)`로 경로를 지정해야 합니다.

```python
from propcalc.mixture import Mixture
from propcalc.backends import get_backend
from propcalc import performance

mix = Mixture.from_wt_pct({"KNO3": 65, "SORBITOL": 35})
res = get_backend("cea").run(mix, Pc_MPa=6.895, Pe_MPa=0.101325)   # 또는 area_ratio=...
print(res.assumptions())              # 평형/동결, Pc, Pe, 백엔드, 배합
print(performance.report(res).summary())
```

## 구현 범위

| 단계 | 내용 | 상태 |
|---|---|---|
| 1 | 화학식 파서, 배합 전처리, 화학량론, 산소수지, 원소수지 검사 | 완료 |
| 2 | 원소분석 + HHV 로 ΔfH° 역산, Channiwala–Parikh 상관식 | 완료 |
| 3 | cea / rocketcea 백엔드, KNSB 65/35 회귀 기준값 | 완료 |
| 4 | c*, Cf, F, Isp, 이상기체 해석해 비교, 경험 효율 옵션 | 완료 |
| 5 | 산화제 스윕, 연료 대체, 혼합 연료, ΔfH° 민감도, 그래프 | 완료 |
| 6 | 교육용 자체 깁스 솔버 (선택) | 미구현 |
| 7 | 내탄도 (선택) | 미구현 |

연소속도 (a, n)는 계산으로 만들 수 없는 실측 물성이므로, 7단계를 구현할 때에도 실측값이 없으면 실행을 거부해야 합니다.

## 데이터 출처 규칙

모든 반응물 값은 `source`와 `verification`(`primary` / `secondary` / `doc_only`)을 가집니다. `primary`가 아닌 입력값은 계산 결과의 `Result.warnings`에 자동으로 실립니다. 현재 DB의 `primary` 항목은 소르비톨 하나뿐입니다.

| 성분 | ΔfH° (kJ/mol) | 검증 수준 | 비고 |
|---|---|---|---|
| KNO3(cr) | −494.63 | secondary | Wagman 1982 원문 미확인. 2차 자료는 −494.5, −494.0 |
| 소르비톨 | −1353.70 ± 1.40 | primary | NIST WebBook 직접 확인(2026-10-30). WebBook이 인용한 Gerasimov 1985 원논문은 미확인 |
| 셀룰로오스(면), (목재), Avicel, 물 | 지시서 기재값 | doc_only | 원문 미확인 |

KNO3의 불확실도는 미확인(null)이라, 민감도 분석에서 쓰려면 `u_kJ_per_mol`을 직접 지정해야 합니다.

바이오매스 연료(코코피트 등)는 현재 DB에 없습니다. 검토한 문헌(코코넛 껍질 분말 논문, 코코피트 퇴비화 논문)은 원소분석이 오류로 보이거나 아예 없어서 채택하지 않았습니다. 이후 후보 문헌의 원문 확인 결과는 [docs/coir_pith_literature_review.md](docs/coir_pith_literature_review.md)에 출처 등급과 함께 있으며, 코코피트의 같은 시료 원소분석·회분·수분·실측 HHV를 원문에서 확인한 문헌은 아직 없습니다(원문을 열람한 3편은 대상 물질이 코코피트가 아니고, 2편은 본문을 열람하지 못했습니다). 새 연료는 `hess.dHf_from_dry_basis_report`로 건조 기준 분석 결과(회분 포함, 실측 HHV)를 받아 입력하며, 의뢰 항목과 보고 형식은 [docs/coir_pith_analysis_request.md](docs/coir_pith_analysis_request.md)에 있습니다.

## 지시서와 다르거나 새로 확인한 점

1. CEA 내장 `KNO3`는 고체가 아닙니다. 이 종의 298 K 엔탈피는 약 −316 kJ/mol로 기체상 값입니다(고체는 −494.6). 이를 그대로 쓰면 KNSB Tc가 1971 K, c*가 1056 m/s로 나옵니다. 지시서의 "KNO3(cr): CEA 내장 DB 사용"은 성립하지 않으므로 고체 KNO3를 커스텀 반응물로 정의합니다. 이 방식에서 Tc 1599 K, c* 910.7 m/s로 문헌 범위에 들어옵니다.
2. `cea.Reactant` 이름은 15자 이하여야 합니다(16자부터 `CEA_INVALID_SIZE`). 백엔드가 `R0_<키>` 형태로 줄여 전달합니다.
3. CEA는 `converged=True`를 보고하면서도 비물리적 값을 낼 수 있습니다. 예: KNO3 87% 이상의 소르비톨(c* 0.1 m/s, M=0), 리그닌(가상 시료) KNO3 30~34%(Isp 3.4×10⁷ s, NaN). `check_physical`이 이를 `ConvergenceError`로 거부하며, 스윕은 해당 행을 지우지 않고 `ok=False`와 사유를 남깁니다.
4. 동결 유동(`mode="frozen"`)은 응축상이 있는 계(KNSB)에서 CEA가 "응축상 온도 범위를 50 K 이상 벗어남"으로 중단하므로 `ConvergenceError`가 납니다. rocketcea도 0.0을 반환합니다. rocketcea 백엔드는 평형 전용입니다.
5. 첫 세션에서는 `webbook.nist.gov`와 `nasa.github.io`에 접속하지 못했으나 2026-10-30 세션에서 둘 다 확인했습니다. `n_frz=2`가 목부터 동결이라는 것은 nasa.github.io/cea의 Example 12("composition frozen at the throat", `n_frz = 2`)와 그 출력(연소실 열만 조성이 다르고 목 이후 열의 질량분율이 같음)으로 확인했습니다. `n_frz=1`이 연소실부터 동결이라는 것은 웹 문서에 명시가 없어 여전히 샘플 주석과 rocketcea 문서 문자열에서 추론한 값입니다. 이 코드는 `frozen_at="chamber"`일 때만 `n_frz=1`을 씁니다.

## 검증 결과 (KNSB 65/35, Pc 6.895 MPa, Pe 0.101325 MPa, 평형)

| 항목 | 이 계산(cea) | rocketcea | 지시서 문헌 범위 |
|---|---|---|---|
| Tc (K) | 1599.37 | 1599.37 | 1550~1650 |
| c* (m/s) | 910.65 | 910.65 | 880~950 |
| 응축상 질량분율 | 0.436 | 0.436 | 0.40~0.47 |
| 유효 분자량 M (g/mol) | 39.87 | 39.87 | 37~42 |
| Isp, Pa=Pe (s) / 진공 (s) | 153.2 / 167.8 | 동일 | 145~170 |

- 두 백엔드는 같은 열역학 데이터에 기반하므로 일치가 독립 검증은 아닙니다. 서로 다른 구현의 교차 확인일 뿐입니다.
- 연소실 조성은 Nakka 평형식(소르비톨 1 mol 당)과 CO2 1.867/1.870, CO 2.489/2.490, H2O 4.827/4.828, H2 2.142/2.145, N2 1.673/1.672, K2CO3 1.644/1.644, KOH 0.053/0.057로 소수 셋째 자리까지 가깝습니다. 지시서가 우려한 "칼륨 화합물 DB 차이로 인한 수 % 차이"는 이 조건에서 관측되지 않았습니다.
- 이상기체 해석해 c*에 CEA의 M_eff와 γ_s를 넣으면 평형 c*와 +0.01%로 일치하지만, CEA가 응축상을 무슬립·열평형 유효 기체로 다루는 것과 같은 가정이므로 독립 검증이 아닙니다. 응축상을 한 분자로 세는 MW(35.41)를 넣으면 +6.1% 과대평가됩니다.

## 5단계 결과 요약

그림과 표는 `notebooks/figures/`에 있습니다(`knsb_sweep`, `fuel_comparison`, `sensitivity`, 각각 PNG와 CSV). 모든 그림 하단에 가정 조건이 표기됩니다.

1. KNSB 산화제 스윕: Isp 최대는 KNO3 72%(155.39 s)로 화학량론(74.27%)보다 연료 과잉 쪽입니다. 72%와 74%를 비교하면 Tc는 −4.0%, M_eff는 −4.6%(47.4 대 49.7)이고 응축상 질량분율은 72%가 0.412, 74%가 0.392로 +2.0%p입니다. √(Tc/M)가 약 +0.3% 늘어 Isp가 +0.15% 높습니다. 최대점은 매우 평평합니다(71%, 72%, 73%에서 155.37, 155.39, 155.31 s). 따라서 "연료 과잉 쪽"은 이 모델 안에서만 성립하는 약한 결론입니다.
2. 연료 대체(KNO3 비율별 최대 Isp): 셀룰로오스 156.50 s @70%(화학량론 74.96%), 소르비톨 155.39 s @72%, 리그닌(가상 시료) 146.6 s @83%(화학량론 83.84%). 소르비톨과 셀룰로오스의 차이 약 1.1 s는 셀룰로오스 ΔfH° ±30 kJ/mol이 만드는 Isp 편차(KNO3 65%에서 ±1.5 s)와 같은 크기이므로, 이 데이터만으로 우열을 확정할 수 없습니다.
3. 흑연 C(gr): 소르비톨 KNSB 스윕 범위(50~85%)에서는 나타나지 않습니다. 더 연료 과잉인 영역에서는 소르비톨 KNO3 30~47%, 셀룰로오스 30~49%, 리그닌(가상 시료) 35~71%에서 나타나며, 리그닌은 KNO3 65%에서도 흑연이 평형 생성물에 포함됩니다.
4. ΔfH° 민감도(KNO3 65%, 상·하한): 셀룰로오스 ±30 kJ/mol은 Tc −29/+28 K, Isp ±1.5 s, 리그닌(가상 시료) ±89 kJ/mol(HHV ±0.5 MJ/kg에서 전파)은 Tc ±31 K, Isp ±3.5 s, 소르비톨 ±1.4 kJ/mol은 Tc ±1.3 K입니다. KNO3의 ΔfH°가 +1 kJ/mol 바뀌면 Tc는 +3.1 K 바뀌는데, KNO3 비중이 크기 때문에 소르비톨 1 kJ/mol(+0.9 K)보다 영향이 큽니다.
5. 모델 아티팩트: K2CO3의 상전이 plateau(약 1173 K) 근처에서는 평형 모델이 γ_s와 c*를 불연속적으로 만듭니다(KNO3 53%에서 γ_s=0.93, 56~58%에서 c* 꺾임). 스윕은 상 공존이 관측된 온도 ±20 K 안의 점에 `flag`를 붙입니다. 리그닌(가상 시료)의 KNO3 54~66% 구간은 전부 이 plateau에 걸려 있어 해당 구간의 결과는 추가로 주의해야 합니다.

리그닌은 지시서의 참고 조성(C10H11.5O2.9)과 가정 HHV 26.8 ± 0.5 MJ/kg로 만든 가상 시료이며, 실제 시료의 원소분석과 HHV로 교체하기 전에는 정량 결론으로 쓰면 안 됩니다.

## 코코피트 대용 분석 (외피 두 시료, 코코피트 아님)

코코피트 문헌을 원문에서 확인하지 못해 DB에는 코코피트가 없습니다. 대신 RMUTP 2025 논문(어린 코코넛 외피, [docs/coir_pith_literature_review.md](docs/coir_pith_literature_review.md))에서 같은 산지 행의 원소분석과 실측 HHV를 읽어 시료별로 `YOUNG_COCONUT_HUSK_RB`, `YOUNG_COCONUT_HUSK_SS`를 등록했습니다. `python notebooks/coir_pith_proxy_analysis.py`로 두 시료의 ΔfH° 한계값 범위에서 성능을 계산합니다.

1. 두 시료의 ΔfH°는 −1463.7 ± 100.4, −1114.7 ± 103.8 kJ/mol로 약 349 kJ/mol 다릅니다. 한 시료 안의 HHV 반복 SD(0.13~0.26 MJ/kg)보다 시료 간 차이(0.87 MJ/kg, 표 기준)가 큽니다.
2. KNO3 65%에서 한계값 전체 범위는 Tc 1286~1687 K(폭 401 K), Isp 139.3~157.5 s입니다. 같은 조건의 소르비톨(153.2 s)이 이 범위 안에 있어 우열을 판단할 수 없습니다.
3. 표의 기준(건조/받은 상태)이 논문에 명시돼 있지 않아 수분 포함 기준으로 가정하고 수분의 H, O를 뺀 뒤 환산했습니다. 표를 건조 기준으로 보면 ΔfH°가 147~196 kJ/mol 바뀝니다(RB −1610.5, SS −1311.0). 이 차이의 반폭(73~98 kJ/mol)을 HHV 반복 SD와 직교 합성해 u에 반영했습니다.
4. 코코피트가 이 외피와 비슷하다는 가정은 검증되지 않았습니다. 위 범위는 "이 두 시료가 보여 주는 규모"일 뿐 코코피트의 범위가 아닙니다.
5. 방침: 코코피트의 같은 시료 원소분석과 실측 HHV를 확보할 방법이 없어, 코코피트는 정량 결론에서 제외합니다. 이 분석은 "코코넛 계열 연료의 ΔfH° 불확실도가 성능에 주는 규모"를 보이는 민감도 사례로만 인용합니다. 보고서에는 코코피트의 Tc, c*, Isp 수치를 쓰지 않습니다.
6. 회분(건조 환산 5~6%)은 반응물에서 제외됩니다. S가 포함된 조성은 화학량론 계산(`stoichiometry`)을 지원하지 않아 이 분석에는 화학량론 산화제 비율이 없습니다.

## 알려진 한계

1. 평형 계산은 무한 반응속도, 단열, 1D 등엔트로피, 입자-기체 완전 평형 가정이므로 Isp의 상한을 줍니다. 학술 문헌상 KNSB 실측 해면 Isp는 이상값의 약 78~90% 수준으로 보고됩니다. 경험 효율 `ETA_ISP_KNSB_AMATEUR`(137/151.5)는 아마추어 자료(Nakka) 기반이므로 명시적으로 선택할 때만 적용합니다.
2. 리그닌은 정의된 분자가 아닙니다. 원소분석과 HHV는 같은 시료에서 얻은 값이어야 하며, 다른 논문 값을 섞으면 ΔfH° 오차가 수백 kJ/mol에 이를 수 있습니다.
3. 크라프트 리그닌은 S, N, 회분(Na 등)을 포함할 수 있어 생성물 종을 추가해야 할 수 있습니다. 현재 `hess.py`는 C, H, O, N, S만 지원합니다.
4. 셀룰로오스와 리그닌은 흡습성이 강하므로 `Mixture.with_moisture`로 수분을 별도 성분(H2O(l))으로 입력해야 합니다.
5. 바이오매스 혼합 추진제의 연소속도는 공개 문헌에서 확인되지 않았습니다. 추력-시간 곡선은 실측 없이는 계산하지 않습니다.
6. PROPEP(제한된 종 DB)과 CEA(NASA Glenn DB)는 칼륨 화합물 데이터 차이로 c*에서 수 % 차이가 날 수 있습니다.
7. HHV를 ΔcH°로 취급하므로 정적 연소열과의 차이(Δn_gas·RT, 리그닌 약 −3.5 kJ/mol)는 무시합니다.
8. 반응물 온도는 298.15 K만 지원합니다.

## 참고문헌

- Gordon, S., McBride, B.J. (1994). NASA RP-1311.
- McBride, B.J., Zehe, M.J., Gordon, S. (2002). NASA/TP-2002-211556.
- NASA CEA: https://github.com/nasa/cea, 문서 https://nasa.github.io/cea/
- NIST Chemistry WebBook, D-sorbitol (CAS 50-70-4): https://webbook.nist.gov/cgi/cbook.cgi?ID=C50704&Units=SI&Mask=2
- Gerasimov, P.A., Blokh, E.L., Gubareva, A.I. (1985). Izv. Vyssh. Uchebn. Zaved., Khim. Khim. Tekhnol. 28, 54-56 (WebBook 인용, 원논문 미확인).
- RocketCEA: https://rocketcea.readthedocs.io/en/latest/new_propellants.html
- Jessup, R.S., Prosen, E.J. (1950). J. Res. NBS 44, 387.
- Colbert, J.C. et al. (1981). J. Res. NBS 86, 655.
- Channiwala, S.A., Parikh, P.P. (2002). Fuel 81, 1051.
- Wagman, D.D. et al. (1982). J. Phys. Chem. Ref. Data 11, Suppl. 2.
- Cox, J.D., Wagman, D.D., Medvedev, V.A. (1989). CODATA Key Values for Thermodynamics.
- Nakka, R. KNSB 열화학 자료 (아마추어 자료, 범위 비교용): https://www.nakka-rocketry.net/sorbchem.html
- J. Chil. Chem. Soc. 68(1), 2023.
