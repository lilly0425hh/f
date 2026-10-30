# 코코피트 문헌 원문 확인 기록 (2026-10-30)

목적: 코코피트(coir pith)의 원소분석(C, H, N, S, O), 회분, 수분, 실측 HHV가 같은 시료로 보고된 문헌이 있는지 원문에서 확인한다. 기준은 [coir_pith_analysis_request.md](coir_pith_analysis_request.md)와 README의 데이터 출처 규칙이다. 확인하지 못한 값은 DB에 넣지 않는다.

결론: 후보 5편 중 같은 시료의 원소분석·회분·수분·실측 HHV를 이 세션에서 원문으로 확인한 것은 0편입니다. 따라서 `ingredients.yaml`에 코코피트를 추가하지 않았습니다.

## 출처 등급

| 등급 | 뜻 | DB 입력 |
|---|---|---|
| A | 이 세션에서 원문 표를 직접 열람했고, 같은 시료의 필요 항목이 모두 있으며 기준(건조/받은 상태)이 명시됨 | 가능 (`primary`) |
| B | 원문을 열람했으나 대상 물질이 다르거나 필요 항목이 빠져 있음 | 불가 |
| C | 검색 결과 목록 등 2차 정보만 봄. 원문 미열람 | 불가 |
| D | 접속 실패. 아무것도 확인하지 못함 | 불가 |

## 후보 문헌별 결과

| 후보 | 등급 | 확인한 내용 | 확인하지 못한 내용 |
|---|---|---|---|
| coirboard.gov.in Coir-Pith.pdf | D | 없음. 프록시 터널이 접속 중 끊겼고(`ws_closed_mid_exchange`), WebFetch는 HTTP 503을 반환 | 문서의 성격, 수치, 인용한 1차 출처 전부 |
| PMC7718462 | B | 아래 상세 | 해당 없음 (코코피트 논문이 아님) |
| chembioagro 10.1186/s40538-023-00431-6 | D | 같은 출판사 도메인으로 301 이동 후 `idp.springer.com` 인증 리다이렉트(303)에서 막힘. `idp.springer.com`은 프록시가 403으로 거부 | 제목, Table 1을 포함한 본문 전부 |
| 태국 RMUTP 코코넛 부산물 논문 | C | 검색 결과 목록에서 "Biofuel and Charcoal Properties of Discarded Aromatic Coconut Residues from Agro-industry"(RMUTP 저널, ph02.tci-thaijo.org 논문 255357)가 후보로 나옴. 목록 설명에 봄브 열량계 HHV와 C, H, O, N, S 분석이 있다고 적혀 있음 | 수치 전부. 이 논문이 코코피트를 포함하는지, 회분·수분이 같은 시료로 있는지. `ph02.tci-thaijo.org`는 이 세션의 네트워크 허용 목록에 없어 열람 불가 |
| Springer 10.1007/s12649-018-0362-2 | D | `idp.springer.com` 인증 리다이렉트(303)에서 막힘 | 제목을 포함해 전부 |

RMUTP 항목은 "태국 RMUTP 코코넛 부산물 논문"이라는 설명과 검색 결과가 일치한다고 제가 추정한 것입니다. 의도하신 논문이 이것인지 확인이 필요합니다.

### PMC7718462 상세 (등급 B)

원문(PMC 전문 HTML)을 열람했습니다. 이 논문은 코코피트가 아닙니다.

- 서지: Cabrales H., Arzola N., Araque O. 외, "The effects of moisture content, fiber length and compaction time on African oil palm empty fruit bunches briquette quality parameters", Heliyon 2020;6(12):e05607, doi 10.1016/j.heliyon.2020.e05607. 시료는 아프리카 기름야자 빈과송이(EFB)입니다. 본문에 coir, coco, pith라는 단어가 없습니다.
- Table 9는 Sheng and Azevedo 상관식 `19.914 − 0.2324·%Ash`로 구한 HHV 추정값 18.607 MJ/kg(상관계수 0.625)입니다. 실측값이 아니라 회분만으로 계산한 추정값입니다.
- Table 6은 C, H, N만 있고(C 40.81/47.69/37.76, H 5.67/6.20/5.97, N 1.01/1.04/1.01 %, 세 산지), 합이 44.75~54.93%라 나머지 약 45~55%는 보고되지 않았습니다. S와 O 값은 표에 없고 본문에 S가 1% 미만이라는 서술만 있습니다.
- Table 7의 실측 HHV는 17.84~18.63 MJ/kg이고, 본문의 회분 5.62%는 세 산지의 평균입니다. 산지별 같은 시료의 회분이 아닙니다.

의뢰하신 후보 목록에 이 논문이 "table 9"와 함께 들어 있는데, 코코피트 자료로는 부적합합니다. 인용 번호가 잘못 옮겨졌을 가능성이 있으니 의도하신 논문의 PMC 번호를 다시 확인해 주세요.

## 열역학 데이터 확인

### 소르비톨 고체 ΔfH° (NIST WebBook, 등급 A)

2026-10-30에 `https://webbook.nist.gov/cgi/cbook.cgi?ID=C50704&Units=SI&Mask=2`를 직접 열람했습니다.

| 항목 | 값 | 방법 | 인용 |
|---|---|---|---|
| ΔfH°(solid) | −1353.7 ± 1.4 kJ/mol | Ccb | Gerasimov, Blokh, et al., 1985 (편집 주석 ALS) |
| ΔcH°(solid) | −3009.4 ± 1.4 kJ/mol | Ccb | 같은 논문. 주석: "Corresponding ΔfH°solid = −1352.5 kJ/mol (simple calculation by NIST; no Washburn corrections)" |
| Cp(solid, 298.15 K) | 241.43 J/(mol·K) | | Lian, Chen, et al., 1982 |

DB의 기존 값 −1353.70 ± 1.40과 일치합니다. 그래서 `SORBITOL`의 `verification`을 `secondary`에서 `primary`로 올리고 출처를 WebBook 직접 확인으로 바꿨습니다. 한계는 두 가지입니다. WebBook이 인용한 Gerasimov 1985 원논문(비타민 C 합성 중간체의 표준 생성열 보고)은 읽지 못했습니다. 또한 같은 페이지의 연소열 환산값(−1352.5)과 1.2 kJ/mol 차이가 있고, 1985년 단일 연소열 측정이어서 ±1.4가 실제 불확실도를 과소평가할 가능성은 배제하지 못합니다. 이 차이는 소르비톨 ±1.4 kJ/mol이 Tc ±1.3 K를 만든다는 README의 민감도 결과 안에서는 작지만, 불확실도를 키워 다시 돌릴지는 판단이 필요합니다.

### RocketSolver `n_frz` (nasa.github.io/cea, 부분 확인)

- 확인됨: Example 12 본문에 "a variable n_frz that we will use later to declare which rocket station should begin the frozen composition"이라 적혀 있고, 같은 예제가 "composition frozen at the throat"이며 `n_frz = 2`를 씁니다. 출력의 질량분율은 첫 열(연소실)만 다르고 둘째 열(목) 이후 10개 열이 모두 같습니다(CO 0.076919 → 0.067456 이후 동일). 따라서 `n_frz=2`는 "목의 평형 조성에서 조성을 고정, 목 이후 동결"이며, 목까지는 평형으로 계산됩니다.
- 확인 못함: `n_frz=1`이 연소실부터 동결이라는 것. 웹 문서의 Python/MATLAB API 설명은 `n_frz: int | None` 시그니처만 있고 값의 의미를 정의하지 않습니다. 이 값은 이전 세션처럼 샘플 주석과 rocketcea 문서 문자열에서 추론한 것으로 남아 있습니다.
- 코드의 `frozen_at="chamber"` 경로(`n_frz=1`)는 이 추론에 의존합니다. `nasa/cea` 저장소 소스는 이 세션의 GitHub 접근 범위 밖이어서 열람하지 않았습니다.

## 다음 단계에 필요한 것

1. 네트워크 허용: `coirboard.gov.in`, `idp.springer.com`(Springer 두 편의 인증 리다이렉트), `ph02.tci-thaijo.org`를 환경 설정의 Network access에서 허용하면 원문을 다시 열람할 수 있습니다. 허용 방법은 https://code.claude.com/docs/en/cloud-environments#network-access 에 있습니다. 허용하지 않고 PDF를 저장소에 올려 주셔도 됩니다.
2. 의도하신 PMC 논문 번호와 RMUTP 논문 제목의 확인.
3. 원문을 열람한 뒤에도 A 등급 조건(같은 시료, 기준 명시, 합 100±1.5%)을 만족하는 문헌이 없을 수 있습니다. 문헌값은 공개 보고서마다 시료의 처리 이력(세척, 완충)이 다르므로, 의뢰서 방식의 직접 분석이 여전히 가장 안전한 경로입니다.
