'use strict';
const assert=require('node:assert/strict'),fs=require('node:fs'),path=require('node:path');
const {pathToFileURL}=require('node:url');
const {chromium}=require(process.env.PLAYWRIGHT_MODULE||'playwright');
const qa=path.join(__dirname,'.qa');fs.mkdirSync(qa,{recursive:true});
const data=JSON.parse(fs.readFileSync(path.join(__dirname,'catalogue.json'),'utf8'));
const records=data.records.filter(r=>r.scope.symmetrizer==='positive_diagonal');
const url=pathToFileURL(path.join(__dirname,'index.html')).href;
(async()=>{
  const browser=await chromium.launch({headless:true,...(process.env.BROWSER_EXECUTABLE?{executablePath:process.env.BROWSER_EXECUTABLE}:{})});
  try{
    const context=await browser.newContext({offline:true,viewport:{width:1440,height:1050},acceptDownloads:true});
    const page=await context.newPage(),errors=[];page.on('pageerror',e=>errors.push(String(e)));
    async function route(id,tab='matrices'){
      await page.goto(url+'#'+id+'/'+tab);
      await page.waitForFunction(([id,tab])=>document.querySelector('#record-header .eyebrow').textContent.endsWith(id)&&document.querySelector('#tab-'+tab).getAttribute('aria-selected')==='true',[id,tab]);
    }
    await route('s2-c01');await page.locator('#navigator-filters summary').click();await page.selectOption('#symmetrizer-filter','positive_diagonal');
    for(const rank of [...new Set(records.map(r=>r.rank))]){
      await page.click(`[data-rank="${rank}"]`);
      await page.waitForFunction(rank=>document.querySelector('#record-header .eyebrow').textContent.startsWith('Rank '+rank),rank);
      assert.equal(await page.locator('#symmetrizer-filter').inputValue(),'positive_diagonal');
      assert.equal(await page.locator('.record-link').count(),records.filter(r=>r.rank===rank).length);
    }
    for(const r of records){
      await route(r.id,'mutation');
      assert.match(await page.locator('#record-header .facts').textContent(),new RegExp('diag \\('+r.datum.symmetrizer.join(', ')+'\\)'));
      await page.locator('[data-vertex="0"]').click();
      assert.match(await page.locator('#graph-caption').textContent(),new RegExp('Vertex 1, weight '+r.slice.symmetrizer[0]+':'));
      const text=await page.locator('#quiver').textContent();
      for(let i=0;i<r.slice.B.length;i++)for(let j=0;j<r.slice.B.length;j++){
        if(r.slice.B[i][j]>0&&r.slice.B[i][j]!==-r.slice.B[j][i])assert(text.includes('('+r.slice.B[i][j]+', '+(-r.slice.B[j][i])+')'),r.id+' valuation');
      }
      await page.click('#tab-exponents');
      assert.match(await page.locator('#panel-exponents .formula').textContent(),/log f = K∨ log\(1 − f\),\s+K∨ = D−1KD/);
    }
    for(const r of records.filter(r=>['certified-finite-orbit','certified-standard','certified-mutation-infinite'].includes(r.notes.quiver.status))){
      await route(r.id,'notes');
      if(r.notes.quiver.certificate.orbit_modulo_relabeling){
        assert((await page.locator('#panel-notes').textContent()).includes(r.notes.quiver.certificate.canonical_orbit_size+' valued exchange matrices up to relabelling'));
        assert(!(await page.locator('#panel-notes').textContent()).includes('undefined'));
      }
      if(r.notes.quiver.geometric_identification)assert((await page.locator('#panel-notes').textContent()).includes('Orbifold type'));
      await page.click('#replay-class');
      await page.waitForFunction(()=>document.querySelector('#mutation-mode').value==='certificate');
      await page.evaluate(()=>{let limit=1000;while(!document.querySelector('#loop-next').disabled&&limit--)document.querySelector('#loop-next').click();if(limit<=0)throw Error('Replay exceeded limit');});
      const actual=await page.locator('#exchange-matrix table').evaluate(t=>[...t.rows].map(row=>[...row.cells].map(c=>Number(c.textContent.replace(/−/g,'-')))));
      assert.deepEqual(actual,r.notes.quiver.certificate.target_B);
    }
    await route('s4-c05','notes');assert(await page.locator('a[href="#s5-c11/notes"]').count());
    await page.screenshot({path:path.join(qa,'weighted-notes.png'),fullPage:true});
    const folded=records.find(r=>r.notes.family.identifications.some(m=>m.category==='Fold'));
    if(folded){
      const note=folded.notes.family.identifications.find(m=>m.category==='Fold');
      await route(folded.id,'notes');
      const detail=page.locator('details.family-note').filter({hasText:note.label}).first();
      await detail.locator('summary').click();
      assert.match(await detail.textContent(),/The folded pair gives this representative/);
      assert(await detail.locator('a[href="#'+note.folding.parent_record+'/notes"]').count());
      await page.screenshot({path:path.join(qa,'weighted-fold-note.png'),fullPage:true});
    }
    await route('s5-c40','matrices');
    const download=page.waitForEvent('download');await page.click('#export-tex');const tex=await download;await tex.saveAs(path.join(qa,'weighted.tex'));
    assert(fs.readFileSync(path.join(qa,'weighted.tex'),'utf8').includes('\\[ D = \\operatorname{diag}('));
    await page.screenshot({path:path.join(qa,'weighted-matrices.png'),fullPage:true});
    await route('s5-c40','family');await page.fill('#scale','2');
    const liftDownload=page.waitForEvent('download');await page.click('#export-lift');const lift=await liftDownload;await lift.saveAs(path.join(qa,'weighted-lift.json'));
    assert.deepEqual(JSON.parse(fs.readFileSync(path.join(qa,'weighted-lift.json'),'utf8')).datum.symmetrizer,data.records.find(r=>r.id==='s5-c40').datum.symmetrizer);
    await route('s5-c40','exponents');await page.screenshot({path:path.join(qa,'weighted-exponents.png'),fullPage:true});
    const largest=records.reduce((a,b)=>a.slice.vertices>b.slice.vertices?a:b);
    await route(largest.id,'exponents');await page.screenshot({path:path.join(qa,'weighted-largest-exponents.png'),fullPage:true});
    const rank6=records.filter(r=>r.rank===6);
    if(rank6.length){
      const maxMultiplicity=r=>Math.max(...r.exponents.multiplicities.map(x=>x.multiplicity));
      const repeated=rank6.reduce((a,b)=>maxMultiplicity(a)>maxMultiplicity(b)?a:b);
      const root=repeated.exponents.multiplicities.find(x=>x.multiplicity===maxMultiplicity(repeated));
      await route(repeated.id,'exponents');await page.locator(`[data-root="${root.m}"]`).click();
      assert((await page.locator('#spectrum-selection').textContent()).includes('multiplicity '+root.multiplicity));
      await page.screenshot({path:path.join(qa,'weighted-rank6-multiplicities.png'),fullPage:true});
      const orbifold=rank6.find(r=>r.notes.quiver.geometric_identification);
      if(orbifold){await route(orbifold.id,'notes');await page.screenshot({path:path.join(qa,'weighted-rank6-orbifold-note.png'),fullPage:true});}
    }
    await page.setViewportSize({width:390,height:844});
    for(const tab of ['matrices','mutation','exponents']){
      await route(largest.id,tab);await page.locator('#panel-'+tab).scrollIntoViewIfNeeded();
      assert(await page.evaluate(()=>document.documentElement.scrollWidth<=window.innerWidth),tab+' mobile width');
      await page.screenshot({path:path.join(qa,'weighted-mobile-'+tab+'.png'),fullPage:true});
    }
    assert.deepEqual(errors,[]);console.log(`PASS: ${records.length} weighted records, rank filters, valuations, dual fixed-point convention, mutation-class replay, crosslinks, exports and mobile layout.`);
  }finally{await browser.close();}
})().catch(e=>{console.error(e);process.exitCode=1;});
