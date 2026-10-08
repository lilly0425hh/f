// database.rules.json 검증: Firebase Realtime Database 에뮬레이터에서 실제 규칙을 실행한다.
// 실행: (에뮬레이터 실행 중) node tests/rules.test.mjs   — README의 "테스트" 절 참고
import fs from 'node:fs';
import { initializeTestEnvironment, assertSucceeds, assertFails } from '@firebase/rules-unit-testing';

const rules = fs.readFileSync(new URL('../database.rules.json', import.meta.url), 'utf8');
const env = await initializeTestEnvironment({ projectId: 'demo-classseat', database: { host: '127.0.0.1', port: 9000, rules } });
const H = 'a'.repeat(64), H2 = 'b'.repeat(64);
const TS = { '.sv': 'timestamp' }, INC = (n) => ({ '.sv': { increment: n } });
let pass = 0, fail = 0;
const t = async (name, p, expectOk = true) => {
  try { await (expectOk ? assertSucceeds(p) : assertFails(p)); pass++; console.log('  ✅', name); }
  catch (e) { fail++; console.log('  ❌', name, '\n     ', String(e.message || e).split('\n')[0]); }
};
const seed = async () => {
  await env.clearDatabase();
  await env.withSecurityRulesDisabled(async (c) => {
    const db = c.database();
    await db.ref('admins/T').set(true);
    await db.ref('roomIndex/r1').set({ name: 'x' });
    await db.ref('pinKeys/r1').set({ [H]: true });
    await db.ref('rooms/r1').set({
      config: { rows: 6, cols: 6, budget: 100, seats: { r1c1: true, r1c2: true, r2c1: true } },
      state: { phase: 'betting', round: 1 },
      roster: { s1: { name: '김가' }, s2: { name: '이나' }, s3: { name: '박다' } },
      constraints: [{ a: 's1', b: 's2' }], audit: { seed: 1 },
    });
  });
};
const ctx = (uid) => env.authenticatedContext(uid).database();
const admin = () => env.authenticatedContext('T').database();
const join = async (uid, sid, name) => { // 정상 입장 절차
  await ctx(uid).ref(`rooms/r1/members/${uid}`).set({ k: H });
  await ctx(uid).ref(`rooms/r1/bindings/${sid}`).set({ uid, name });
  await ctx(uid).ref(`rooms/r1/uids/${uid}`).set(sid);
};
const bet = (seats) => ({ seats, confirmed: true, submittedAt: TS });

console.log('\n[1] PIN 검증과 방 입장');
await seed();
await t('목록 조회 불가: pinKeys/r1 읽기 실패', ctx('u1').ref('pinKeys/r1').get(), false);
await t('정확한 키 읽기 성공', ctx('u1').ref(`pinKeys/r1/${H}`).get());
await t('pinKeys 전체 읽기 실패', ctx('u1').ref('pinKeys').get(), false);
await t('PIN 검증 전: config 읽기 실패', ctx('u1').ref('rooms/r1/config').get(), false);
await t('틀린 PIN 해시로 members 쓰기 실패', ctx('u1').ref('rooms/r1/members/u1').set({ k: H2 }), false);
await t('형식이 틀린 k 쓰기 실패', ctx('u1').ref('rooms/r1/members/u1').set({ k: 'zz' }), false);
await t('남의 uid로 members 쓰기 실패', ctx('u1').ref('rooms/r1/members/u2').set({ k: H }), false);
await t('올바른 PIN 해시로 members 쓰기 성공', ctx('u1').ref('rooms/r1/members/u1').set({ k: H }));
await t('member 후 config 읽기 성공', ctx('u1').ref('rooms/r1/config').get());
await t('member 후 roster/state/totals 읽기 성공', Promise.all(['roster', 'state', 'totals'].map((p) => ctx('u1').ref(`rooms/r1/${p}`).get())));
await t('다른 방(r2) 읽기 실패', ctx('u1').ref('rooms/r2/config').get(), false);
await t('constraints 읽기 실패(교사 전용)', ctx('u1').ref('rooms/r1/constraints').get(), false);
await t('audit 읽기 실패(교사 전용)', ctx('u1').ref('rooms/r1/audit').get(), false);
await t('rooms/r1 전체 읽기 실패', ctx('u1').ref('rooms/r1').get(), false);
await t('rooms 전체 읽기 실패', ctx('u1').ref('rooms').get(), false);
await t('config 쓰기 실패', ctx('u1').ref('rooms/r1/config/budget').set(9999), false);
await t('roster 쓰기 실패', ctx('u1').ref('rooms/r1/roster/s9').set({ name: 'x' }), false);
await t('state.phase 쓰기 실패', ctx('u1').ref('rooms/r1/state/phase').set('done'), false);
await t('roomIndex 읽기 실패', ctx('u1').ref('roomIndex').get(), false);
await t('다른 사람 admins 읽기 실패', ctx('u1').ref('admins/T').get(), false);
await t('스스로 admins 등록 실패', ctx('u1').ref('admins/u1').set(true), false);

console.log('\n[2] 학번 binding');
await t('명단에 없는 학번 binding 실패', ctx('u1').ref('rooms/r1/bindings/s9').set({ uid: 'u1', name: 'x' }), false);
await t('이름이 다르면 binding 실패', ctx('u1').ref('rooms/r1/bindings/s1').set({ uid: 'u1', name: '틀림' }), false);
await t('남의 uid로 binding 실패', ctx('u1').ref('rooms/r1/bindings/s1').set({ uid: 'u2', name: '김가' }), false);
await t('정상 binding 성공', ctx('u1').ref('rooms/r1/bindings/s1').set({ uid: 'u1', name: '김가' }));
await t('uids 정상 기록 성공', ctx('u1').ref('rooms/r1/uids/u1').set('s1'));
await t('내 binding 읽기 성공', ctx('u1').ref('rooms/r1/bindings/s1').get());
await ctx('u2').ref('rooms/r1/members/u2').set({ k: H });
await t('다른 uid가 같은 학번 binding 실패', ctx('u2').ref('rooms/r1/bindings/s1').set({ uid: 'u2', name: '김가' }), false);
await t('남의 binding 읽기 실패(이미 묶인 학번)', ctx('u2').ref('rooms/r1/bindings/s1').get(), false);
await t('비어 있는 학번 binding 읽기 허용(null)', ctx('u2').ref('rooms/r1/bindings/s2').get());
await t('binding 목록 읽기 실패', ctx('u2').ref('rooms/r1/bindings').get(), false);
await t('binding 삭제(초기화)는 학생 불가', ctx('u1').ref('rooms/r1/bindings/s1').remove(), false);
await t('같은 uid가 두 번째 학번 binding 실패', ctx('u1').ref('rooms/r1/bindings/s2').set({ uid: 'u1', name: '이나' }), false);
await t('남의 학번을 내 uids에 기록 실패', ctx('u2').ref('rooms/r1/uids/u2').set('s1'), false);
await t('다른 사람 uids 읽기 실패', ctx('u2').ref('rooms/r1/uids/u1').get(), false);

console.log('\n[3] 베팅(bets) 쓰기·읽기');
await t('내 학번에 정상 베팅 성공', ctx('u1').ref('rooms/r1/bets/s1').set(bet({ r1c1: 30, r1c2: 20 })));
await t('빈 베팅(0개 제출) 성공', ctx('u1').ref('rooms/r1/bets/s1').set({ confirmed: true, submittedAt: TS }));
await t('남의 학번에 베팅 실패', ctx('u1').ref('rooms/r1/bets/s2').set(bet({ r1c1: 5 })), false);
await t('음수 토큰 실패', ctx('u1').ref('rooms/r1/bets/s1').set(bet({ r1c1: -5 })), false);
await t('0 토큰 항목 실패(0은 저장하지 않음)', ctx('u1').ref('rooms/r1/bets/s1').set(bet({ r1c1: 0 })), false);
await t('소수 토큰 실패', ctx('u1').ref('rooms/r1/bets/s1').set(bet({ r1c1: 1.5 })), false);
await t('budget 초과 단일 값(101) 실패', ctx('u1').ref('rooms/r1/bets/s1').set(bet({ r1c1: 101 })), false);
await t('비활성/없는 자리 실패', ctx('u1').ref('rooms/r1/bets/s1').set(bet({ r9c9: 5 })), false);
await t('문자열 토큰 실패', ctx('u1').ref('rooms/r1/bets/s1').set(bet({ r1c1: '5' })), false);
await t('정의되지 않은 필드 실패', ctx('u1').ref('rooms/r1/bets/s1').set({ ...bet({ r1c1: 5 }), hack: 1 }), false);
await t('미래 submittedAt 실패', ctx('u1').ref('rooms/r1/bets/s1').set({ seats: { r1c1: 5 }, confirmed: true, submittedAt: Date.now() + 1e9 }), false);
await t('confirmed 가 boolean 아니면 실패', ctx('u1').ref('rooms/r1/bets/s1').set({ seats: { r1c1: 5 }, confirmed: 'yes', submittedAt: TS }), false);
await t('내 베팅 읽기 성공', ctx('u1').ref('rooms/r1/bets/s1').get());
await t('남의 베팅 읽기 실패', ctx('u2').ref('rooms/r1/bets/s1').get(), false);
await t('bets 전체 읽기 실패(학생)', ctx('u1').ref('rooms/r1/bets').get(), false);

console.log('\n[4] 합계(totals) 증감');
await seed(); await join('u1', 's1', '김가'); await join('u2', 's2', '이나');
const root = (uid) => ctx(uid).ref();
await t('베팅+합계 증가 원자적 update 성공', root('u1').update({ 'rooms/r1/bets/s1': bet({ r1c1: 30 }), 'rooms/r1/totals/r1c1': INC(30) }));
await t('두 번째 학생도 같은 자리에 증가 성공', root('u2').update({ 'rooms/r1/bets/s2': bet({ r1c1: 10 }), 'rooms/r1/totals/r1c1': INC(10) }));
const tot = (await admin().ref('rooms/r1/totals/r1c1').get()).val();
await t(`합계가 40(서버 increment 반영): ${tot}`, tot === 40 ? Promise.resolve() : Promise.reject(new Error('합계 ' + tot)));
await t('한 번에 +budget 초과 증가(+101) 실패', root('u1').update({ 'rooms/r1/totals/r1c2': INC(101) }), false);
await t('한 번에 -budget 초과 감소 실패', root('u1').update({ 'rooms/r1/totals/r1c1': INC(-101) }), false);
await t('음수가 되는 감소 실패', root('u1').update({ 'rooms/r1/totals/r1c2': INC(-5) }), false);
await t('문자열 합계 실패', ctx('u1').ref('rooms/r1/totals/r1c2').set('x'), false);
await t('소수 합계 실패', ctx('u1').ref('rooms/r1/totals/r1c2').set(1.5), false);
await ctx('u3').ref('rooms/r1/members/u3').set({ k: H });
await t('입장(binding)하지 않은 member의 합계 쓰기 실패', root('u3').update({ 'rooms/r1/totals/r1c1': INC(5) }), false);
await t('member가 아닌 uid의 합계 쓰기 실패', root('u9').update({ 'rooms/r1/totals/r1c1': INC(5) }), false);
await t('합계 읽기(member) 성공', ctx('u1').ref('rooms/r1/totals').get());
await t('합계 읽기(비member) 실패', ctx('u9').ref('rooms/r1/totals').get(), false);

console.log('\n[5] 단계(phase)별 잠금과 결과 공개');
await admin().ref('rooms/r1/state/phase').set('locked');
await t('locked 이후 베팅 쓰기 실패', ctx('u1').ref('rooms/r1/bets/s1').set(bet({ r1c1: 5 })), false);
await t('locked 이후 합계 쓰기 실패', root('u1').update({ 'rooms/r1/totals/r1c1': INC(1) }), false);
await t('locked 이후 내 베팅 읽기는 가능', ctx('u1').ref('rooms/r1/bets/s1').get());
await admin().ref('rooms/r1/state/phase').set('drawing');
await admin().ref('rooms/r1/result').set({ assign: { s1: 'r1c1' }, seed: 7 });
await t('drawing 중 result 읽기 실패(연출 전 스포일러 방지)', ctx('u1').ref('rooms/r1/result').get(), false);
await admin().ref('rooms/r1/state/phase').set('done');
await t('done 이후 result 읽기 성공', ctx('u1').ref('rooms/r1/result').get());
await t('비member는 done이어도 result 읽기 실패', ctx('u9').ref('rooms/r1/result').get(), false);
await t('학생의 result 쓰기 실패', ctx('u1').ref('rooms/r1/result/assign/s1').set('r1c2'), false);
await t('잘못된 phase 값은 교사도 쓰기 실패', admin().ref('rooms/r1/state/phase').set('hacked'), false);

console.log('\n[6] 교사(admins)');
await seed();
await t('교사: 방 전체 읽기', admin().ref('rooms/r1').get());
await t('교사: roomIndex 읽기', admin().ref('roomIndex').get());
await t('교사: 방 생성(roomIndex/pinKeys/rooms)', admin().ref().update({ 'roomIndex/r2': { name: 'n' }, [`pinKeys/r2/${H2}`]: true, 'rooms/r2/state': { phase: 'setup', round: 1 } }));
await t('교사: 잘못된 형식의 pinKey 실패', admin().ref('pinKeys/r2/not-a-hash').set(true), false);
await t('교사: totals를 한 번에 덮어쓰기(재계산)', admin().ref('rooms/r1/totals').set({ r1c1: 500, r1c2: 3 }));
await t('교사: 방 삭제', admin().ref().update({ 'roomIndex/r2': null, 'pinKeys/r2': null, 'rooms/r2': null }));
await t('교사: bets 전체 읽기', admin().ref('rooms/r1/bets').get());
await t('교사: 학생 binding 초기화', admin().ref('rooms/r1/bindings/s1').remove());
await t('교사: 자신의 admins 읽기', admin().ref('admins/T').get());
await t('일반 사용자: roomIndex 쓰기 실패', ctx('u1').ref('roomIndex/rx').set({ name: 'x' }), false);
await t('일반 사용자: pinKeys 쓰기 실패', ctx('u1').ref(`pinKeys/r1/${H2}`).set(true), false);
await t('비로그인: 모든 읽기 실패', env.unauthenticatedContext().database().ref(`pinKeys/r1/${H}`).get(), false);

console.log('\n[7] 지난 회차 자리(history)');
await seed(); await join('u1', 's1', '김가');
await admin().ref('rooms/r1/history/1').set({ round: 1, layout: { rows: 6, cols: 6, groups: 3 }, assign: { s1: 'r1c1' }, types: { s1: 'bet' } });
await t('member는 history 읽기 성공(어느 단계에서든)', ctx('u1').ref('rooms/r1/history').get());
await t('member는 특정 회차 history 읽기 성공', ctx('u1').ref('rooms/r1/history/1').get());
await t('비member는 history 읽기 실패', ctx('u9').ref('rooms/r1/history').get(), false);
await t('학생의 history 쓰기 실패', ctx('u1').ref('rooms/r1/history/2').set({ round: 2, assign: { s1: 'r1c2' } }), false);
await t('학생의 history 수정 실패', ctx('u1').ref('rooms/r1/history/1/assign/s1').set('r1c2'), false);
await t('학생의 history 삭제 실패', ctx('u1').ref('rooms/r1/history').remove(), false);
await t('교사: 새 회차 일괄 갱신(삭제 + history 보관)', admin().ref().update({ 'rooms/r1/bets': null, 'rooms/r1/result': null, 'rooms/r1/history/2': { round: 2, assign: { s1: 'r1c2' } } }));

console.log(`\n결과: 통과 ${pass} / 실패 ${fail}`);
await env.cleanup();
process.exit(fail ? 1 : 0);
