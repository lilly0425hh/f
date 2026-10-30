# 코코피트 문헌 원문 확인 기록 (2026-10-30, 2차 확인 반영)

목적: 코코피트(coir pith)의 원소분석(C, H, N, S, O), 회분, 수분, 실측 HHV가 같은 시료로 보고된 문헌이 있는지 원문에서 확인한다. 기준은 [coir_pith_analysis_request.md](coir_pith_analysis_request.md)와 README의 데이터 출처 규칙이다. 확인하지 못한 값은 DB에 넣지 않는다.

결론: 후보 5편 중 코코피트의 같은 시료 원소분석·회분·수분·실측 HHV를 원문으로 확인한 것은 0편입니다. 1차 확인에서는 네트워크 제한으로 4편을 열람하지 못했고, 2차 확인에서 3편을 열람했습니다. 열람한 문헌은 모두 대상 물질이 코코피트가 아니었고, coirboard 1편은 여전히 열람하지 못했습니다. 따라서 `ingredients.yaml`에 코코피트를 추가하지 않았습니다.

## 출처 등급

| 등급 | 뜻 | DB 입력 |
|---|---|---|
| A | 원문 표를 직접 열람했고, 코코피트의 필요 항목이 같은 시료로 모두 있으며 기준(건조/받은 상태)이 명시됨 | 가능 (`primary`) |
| B | 원문을 열람했으나 대상 물질이 코코피트가 아니거나 필요 항목이 빠져 있음 | 불가 |
| C | 초록이나 검색 결과 등 일부만 봄. 표를 열람하지 못함 | 불가 |
| D | 접속 실패. 아무것도 확인하지 못함 | 불가 |

## 후보 문헌별 결과

| 후보 | 등급 | 실제 내용 |
|---|---|---|
| coirboard.gov.in Coir-Pith.pdf | D | 호스트는 허용됐으나 서버 쪽에서 연결이 매번 끊김(`ws_closed_mid_exchange`)이고 WebFetch는 503. 문서 내용은 전혀 모름 |
| PMC7718462 | B | Heliyon 2020, 아프리카 기름야자 빈과송이(EFB) 브리켓 논문 |
| chembioagro 10.1186/s40538-023-00431-6 | B | 아몬드 껍질과 플라스틱(PE, PS)의 촉매 열분해 논문 |
| 태국 RMUTP 논문 | B | 방향성 어린 코코넛의 껍질·외피·외피+껍질·빈과송이 분석 논문. 코코피트 아님 |
| Springer 10.1007/s12649-018-0362-2 | C | 코이어 섬유(coir fibers) 펠릿 논문. 유료 구독이라 초록만 열람 |

의뢰하신 후보 설명과 실제 논문이 맞지 않는 곳이 여럿입니다. PMC7718462와 chembioagro는 번호가 잘못 옮겨졌을 가능성이 큽니다. 의도하신 논문이 따로 있다면 PMC 번호나 DOI를 다시 알려 주세요.

### PMC7718462 (등급 B)

Cabrales H., Arzola N., Araque O. 외, "The effects of moisture content, fiber length and compaction time on African oil palm empty fruit bunches briquette quality parameters", Heliyon 2020;6(12):e05607, doi 10.1016/j.heliyon.2020.e05607. 본문에 coir, coco, pith라는 단어가 없습니다.

- Table 9는 Sheng and Azevedo 상관식 `19.914 − 0.2324·%Ash`로 구한 HHV 추정값 18.607 MJ/kg(상관계수 0.625)입니다. 실측값이 아닙니다.
- Table 6은 C, H, N만 있고(C 40.81/47.69/37.76, H 5.67/6.20/5.97, N 1.01/1.04/1.01 %, 세 산지) 합이 44.75~54.93%입니다. S와 O 값은 표에 없습니다.
- Table 7의 실측 HHV는 17.84~18.63 MJ/kg이고, 본문의 회분 5.62%는 세 산지의 평균입니다.

### chembioagro 10.1186/s40538-023-00431-6 (등급 B)

Alcazar-Ruiz A., Sanchez-Silva L., Dorado F., "Hydrocarbon selectivity enhancement through catalytic fast co-pyrolysis of almond shell and plastic blends", Chemical and Biological Technologies in Agriculture, 2023. Springer 사이트의 본문 HTML을 열람했고 coir pith, cocopeat이라는 표현이 없습니다. Table 1은 아몬드 껍질과 PE, PS의 분석표이며, 본문에 "HHV shown in Table 1 was calculated from the elemental analysis ... with the following empirical correlation (Eq. 1)"이라 적혀 있어 HHV가 실측이 아니라 상관식 계산값입니다. 표 본문은 HTML에 렌더링되지 않아 수치는 확인하지 못했습니다. 코코피트가 아니므로 수치를 추가로 확인할 필요는 없습니다.

### Springer 10.1007/s12649-018-0362-2 (등급 C)

Stelte W., Barsberg S.T., Clemons C., Morais J.P.S., de Freitas Rosa M., Sanadi A.R., "Coir Fibers as Valuable Raw Material for Biofuel Pellet Production", Waste and Biomass Valorization, 2019/11. 구독 문서라 초록만 열람했습니다. 초록은 코이어 섬유의 리그닌 약 41~42 wt%, 회분 15 wt% 미만, 발열량 약 18.5~19 MJ/kg(초록 표기 "high heating value")을 적고 있습니다. 원소분석은 초록에 없고, 대상이 섬유(갈색·흰색)라 코코피트와 같은 물질인지도 확인하지 못했습니다. 본문 표를 읽을 수 있을 때까지 어떤 값도 채택하지 않습니다.

### 태국 RMUTP 논문 (등급 B, 원문 PDF 열람)

Pratinthong N., Vinitnantharat S., Sukkhee N., Intarapong P., Manakij A., "Biofuel and Charcoal Properties of Discarded Aromatic Coconut Residues from Agro-industry", RMUTP Research Journal Sciences and Technology 19(1), 181-194, 2025, doi 10.14456/jrmutp.2025.16. 저자는 KMUTT 소속이고 게재지가 RMUTP 저널입니다.

이 논문의 시료는 방향성 어린 코코넛(aromatic coconut)의 외피(เปลือกมะพร้าว, husk), 껍질(กะลา, shell), 빈과송이(ทะลายเปล่า)이며 코코피트(ขุยมะพร้าว)는 없습니다. 형식은 본 작업에 가장 가깝습니다. Table 1(근사분석, HHV)과 Table 2(원소분석)가 같은 산지·시료 행으로 보고되고, 원소분석은 3회 반복 평균±표준편차이며, HHV는 LECO AC-350 봄브 열량계로 측정했습니다. 1차 확인에서 후보로 올렸던 "원소분석과 실측 HHV가 같은 시료로 보고된" 조건은 이 논문이 충족합니다. 다만 다음 한계가 있습니다.

1. 대상 물질: 어린 코코넛 외피이며 코코피트가 아닙니다.
2. 기준 미명시: 표에 "ร้อยละโดยน้ำหนัก"(wt%)라고만 되어 있고 건조/받은 상태를 적지 않았습니다. Table 1은 수분+회분+휘발분+고정탄소=100.00%이므로 수분을 포함한 기준(풍건 상태로 추정)입니다. HHV는 본문에 "건조된 시료"(เศษ...แห้งก่อนเผา)로 적혀 있어 수분 포함 여부를 확정할 수 없습니다. `dHf_from_dry_basis_report`는 건조 기준을 요구하므로 그대로 넣을 수 없습니다.
3. O는 차이법: 식 (2)는 `O = 100 − (C+H+N+S)`로 인쇄돼 있으나, 표의 O는 6행 모두 `100 − (C+H+N+S) − 회분`과 0.01 이내로 일치합니다(예: 외피 Ratchaburi 45.58+6.13+0.34+0.07+41.83+회분 6.05 = 100.00). 인쇄된 식과 표의 계산이 다릅니다. 또한 합이 구성상 항상 100.00이 되므로 의뢰서의 "합 100±1.5%" 점검이 아무것도 걸러내지 못합니다.
4. 같은 시료 여부: 행 이름이 같을 뿐, 같은 균질화 배치인지는 논문에 명시돼 있지 않습니다.

표의 값(외피, Ratchaburi/Samut Songkhram): C 45.58/45.93, H 6.13/5.89, N 0.34/0.36, S 0.07/0.04, O(차이) 41.83/42.88, 회분 6.05/4.90, 수분 3.51/4.73, HHV 17.52/18.39 MJ/kg. 이 값은 코코피트가 아니므로 DB에 넣지 않았습니다. 어린 코코넛 외피를 별도 연료로 다루고 싶으시면 기준 문제(2번)부터 해결해야 하며, 저자에게 기준을 문의하는 것이 가장 정확합니다.

## 열역학 데이터 확인 (1차 확인에서 완료)

### 소르비톨 고체 ΔfH° (NIST WebBook, 등급 A)

`https://webbook.nist.gov/cgi/cbook.cgi?ID=C50704&Units=SI&Mask=2`를 직접 열람했습니다.

| 항목 | 값 | 방법 | 인용 |
|---|---|---|---|
| ΔfH°(solid) | −1353.7 ± 1.4 kJ/mol | Ccb | Gerasimov, Blokh, et al., 1985 (편집 주석 ALS) |
| ΔcH°(solid) | −3009.4 ± 1.4 kJ/mol | Ccb | 같은 논문. 주석: "Corresponding ΔfH°solid = −1352.5 kJ/mol (simple calculation by NIST; no Washburn corrections)" |
| Cp(solid, 298.15 K) | 241.43 J/(mol·K) | | Lian, Chen, et al., 1982 |

DB의 기존 값 −1353.70 ± 1.40과 일치해 `SORBITOL`을 `secondary`에서 `primary`로 올렸습니다. WebBook이 인용한 Gerasimov 1985 원논문은 읽지 못했고, 연소열 환산값(−1352.5)과 1.2 kJ/mol 차이가 있으며, 1985년 단일 측정이어서 ±1.4가 실제 불확실도를 과소평가할 가능성은 배제하지 못합니다.

### RocketSolver `n_frz` (nasa.github.io/cea, 부분 확인)

- 확인됨: Example 12 본문이 "a variable n_frz that we will use later to declare which rocket station should begin the frozen composition"이라 적고, 같은 예제가 "composition frozen at the throat"에 `n_frz = 2`를 씁니다. 출력의 질량분율은 첫 열(연소실)만 다르고 둘째 열(목) 이후 10개 열이 모두 같습니다. 따라서 `n_frz=2`는 목의 평형 조성에서 조성을 고정하는 것입니다.
- 확인 못함: `n_frz=1`이 연소실부터 동결이라는 것. 웹 문서의 API 설명은 `n_frz: int | None` 시그니처만 있고 값의 의미를 정의하지 않습니다.

## 남은 일

1. coirboard Coir-Pith.pdf: 이 세션에서는 서버가 연결을 끊습니다. PDF를 저장소에 올려 주시면 열람하겠습니다.
2. Springer 12649 본문: 구독이 필요합니다. 기관 접근으로 받은 PDF를 올려 주시면 Table의 원소분석과 기준을 확인하겠습니다. 다만 대상이 코이어 섬유이므로 코코피트와 같은 물질이라고 가정하지 않겠습니다.
3. 문헌에서 코코피트의 A 등급 자료가 나오지 않으면, 의뢰서 방식의 직접 분석이 가장 안전한 경로입니다.
