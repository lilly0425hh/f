# 코코피트 문헌 수집 안내 (코코피트 일반 사례 계산용)

목적: 특정 로트가 아니라 "코코피트 일반"에 대한 사례 계산이므로, 한 논문의 한 시료 값을 출처와 함께 쓰는 것으로 충분합니다. 다만 같은 논문의 두 시료도 HHV가 0.87 MJ/kg 다르므로, 결과는 "이 논문의 이 시료"의 값이며 코코피트의 대표값이 아닙니다. 가능하면 서로 다른 논문 2편 이상을 확보해 범위로 보고하세요.

## 채택 조건 (모두 만족해야 A 등급)

1. 시료가 코코피트(coir pith, coir dust, coco peat, coconut pith)이며 섬유(fiber)나 껍질(shell, husk)이 아님.
2. 한 표(한 시료 행)에 C, H, N, S, O가 있음. O가 차이법이면 그렇게 적혀 있음.
3. 같은 시료의 회분, 수분, 봄브 열량계로 측정한 HHV(또는 GCV, 총발열량)가 있음.
4. 기준(건조 / 받은 상태 / 건조 무회분)이 표나 방법 절에 명시됨.
5. 발열량이 상관식으로 계산한 값이 아님.

하나라도 빠지면 쓰지 않습니다. 원소분석은 이 논문, 발열량은 저 논문에서 가져오는 조합은 금지입니다.

## 검색어

영문(Google Scholar, ScienceDirect, Springer, MDPI 등 소속 기관 접속으로): "coir pith" 또는 "coir dust" 또는 "coco peat" 또는 "coconut pith"에 아래 중 하나를 붙입니다.

- pyrolysis, gasification, torrefaction, hydrothermal carbonization, hydrochar, biochar, combustion
- briquette, pellet, "fuel properties", "feedstock characterization", "ultimate analysis", "proximate analysis"

열화학 변환 논문은 원료 특성 표를 앞부분에 두는 경우가 많습니다. 초록이 아니라 본문의 "Materials and methods"와 첫 번째 표를 보세요.

한국어(DBpia, KCI, RISS): "코코피트", "코이어 더스트", "코코넛 피스", "코코넛 부산물 연료", "열분해", "가스화", "고형연료".

## 우선 확인할 후보 (모두 본문 미확인)

1. Air gasification characteristics of coir pith in a circulating fluidized bed gasifier (ScienceDirect, S0973082609000519)
2. Choudhury 외, "Production of bio-oil from coir pith via pyrolysis..." (Biomass Conversion and Biorefinery 2021, doi 10.1007/s13399-020-00630-3)
3. IIT Hyderabad 코코넛 코이어와 피스 열수액화 논문 (Energy 2021, 제목 미확인)
4. coirboard.gov.in/wp-content/uploads/2016/07/Coir-Pith.pdf (검색 요약상 일반 서술일 가능성이 큼)
5. Phyllis2(ecn.nl/phyllis2)에서 "coir", "coconut coir", "coconut husk" 검색. 레코드가 인용한 원문 출처와 발열량의 실측 여부를 확인. DB 자체는 2차 모음이므로 원문 확인 전에는 secondary입니다.

## 가져오실 것

PDF 또는 아래를 채팅에 붙여 주세요.

- 표 전체(캡션 포함)와 그 표의 각주
- 기준(건조/받은 상태)을 설명한 문장
- 산소를 어떻게 구했는지 설명한 문장
- 발열량 측정 방법(봄브 열량계, 표준 번호) 문장
- 서지 정보(저자, 연도, 저널, DOI)

## 받은 뒤 제가 하는 일

1. 위 5개 조건을 대조해 A/B/C/D 등급을 매깁니다.
2. 건조 기준이면 `dHf_from_dry_basis_report`에 그대로 넣습니다. 받은 상태 기준이면 수분의 H, O를 제거하고 건조 환산합니다. 전체를 (1−수분)으로 나누는 환산은 O가 수분을 포함하면 합이 100%를 넘어 거부됩니다(`tests/test_coconut_husk.py`).
3. 값이 모두 원문 표에서 직접 확인되었을 때만 `data/ingredients.yaml`에 시료별로 넣고, 기준 가정과 불확실도의 성격을 `note`에 적습니다.
4. 코코피트 정량 결론은 A 등급 자료가 생겼을 때 README의 제외 방침을 재검토한 뒤 씁니다.
