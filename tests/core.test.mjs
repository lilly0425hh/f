// index.html 안의 "순수 함수 영역"(<draw-core> 표식 사이)만 꺼내 Node에서 자체 검증을 실행한다.
// 실행: node tests/core.test.mjs
import fs from 'node:fs';
const src = fs.readFileSync(new URL('../index.html', import.meta.url), 'utf8');
const m = src.match(/\/\* <draw-core> \*\/([\s\S]*?)\/\* <\/draw-core> \*\//);
if (!m) throw new Error('draw-core 영역을 찾지 못했습니다');
const { selfTest, runDraw, makeRng, areAdjacent, activeSeatIds } = new Function(m[1] + '\nreturn { selfTest, runDraw, makeRng, areAdjacent, activeSeatIds };')();

const out = selfTest({ runs: 1000, log: console.log });

// 추가: 근접 금지를 만족할 수 없는 경우(전원이 서로 금지) → 대체 탐색 후 위반 수를 정직하게 보고하는지
const cfg = { rows: 2, cols: 4, groups: 1, mode: 'pair', disabled: {} };
const students = ['1', '2', '3', '4', '5', '6'];
const forbidden = [];
for (let i = 0; i < students.length; i++) for (let j = i + 1; j < students.length; j++) forbidden.push({ a: students[i], b: students[j] });
const r = runDraw({ students, seats: activeSeatIds(cfg), validBets: {}, forbidden, cfg, seed: 1 });
const okInfeasible = r.stats.fallback === true && r.hardViolations > 0 && Object.keys(r.assign).length === 6;
console.log(`${okInfeasible ? '✅' : '❌'} 만족 불가 제약: 대체 배치 후 위반 ${r.hardViolations}쌍을 보고 (전원 배정 ${Object.keys(r.assign).length}명)`);

// 분단 사이 통로는 인접이 아니다 / 같은 분단의 대각선은 인접이다
const c2 = { rows: 3, cols: 4, groups: 2 };
const okAdj = areAdjacent('r1c2', 'r2c3', c2) === false && areAdjacent('r1c1', 'r2c2', c2) === true && areAdjacent('r1c1', 'r3c1', c2) === false;
console.log(`${okAdj ? '✅' : '❌'} 인접 정의(체비쇼프 거리 1, 분단 통로 제외)`);

process.exit(out.ok && okInfeasible && okAdj ? 0 : 1);
