'use strict';
const assert=require('node:assert/strict'),path=require('node:path'),fs=require('node:fs');
const {pathToFileURL}=require('node:url');
const {chromium}=require(process.env.PLAYWRIGHT_MODULE||'playwright');
(async()=>{
  const browser=await chromium.launch({headless:true,executablePath:process.env.BROWSER_EXECUTABLE});
  try {
    const page=await browser.newPage({viewport:{width:1450,height:1060}}),errors=[];
    page.on('pageerror',e=>errors.push(String(e)));
    const root=pathToFileURL(path.join(__dirname,'index.html')).href;
    const records=JSON.parse(fs.readFileSync(path.join(__dirname,'catalogue.json'),'utf8')).records;
    await page.goto(root+'#r6-c88/matrices');
    assert.equal(await page.locator('.record-link').count(),records.filter(r=>r.rank===6).length);
    await page.screenshot({path:path.join(__dirname,'.qa/rank6-matrices.png'),fullPage:true});
    await page.goto(root+'#r6-c88/exponents');
    await page.screenshot({path:path.join(__dirname,'.qa/rank6-exponents.png'),fullPage:true});
    await page.click('#tab-evidence');
    const links=await page.locator('#panel-evidence a').evaluateAll(xs=>xs.map(x=>x.getAttribute('href')));
    assert(!links.some(x=>x.includes('null')||x.includes('undefined')));
    assert(links.some(x=>/smt_queries(?:-\d+)?\.zip$/.test(x)));
    await page.click('#tab-proof');
    const proof=await page.locator('#panel-proof').textContent();
    assert.match(proof,/108 indecomposable families/);assert.match(proof,/2,887,440/);
    await page.screenshot({path:path.join(__dirname,'.qa/rank6-classification.png'),fullPage:true});
    await page.goto(root+'#r6-c102/notes');
    assert.match(await page.locator('#panel-notes').textContent(),/Sphere with 8 punctures/);
    await page.screenshot({path:path.join(__dirname,'.qa/rank6-surface-note.png'),fullPage:true});
    await page.goto(root+'#r6-c108/exponents');
    assert.match(await page.locator('#record-header').textContent(),/Class 108/);
    for(const id of ['r6-c30','r6-c07','r6-c102']){
      const certificate=records.find(r=>r.id===id).notes.quiver.certificate;
      await page.goto(root+'#'+id+'/notes');await page.click('#replay-class');
      await page.evaluate(length=>{
        for(let i=0;i<length;i++){
          const next=document.querySelector('#loop-next');
          if(next.disabled)throw Error('Certificate stopped before its final mutation');
          next.click();
        }
      },certificate.mutation_path.length);
      assert.match(await page.locator('#step-label').textContent(),/Certificate target reached/);
      const target=await page.locator('#exchange-matrix table').evaluate(t=>[...t.rows].map(row=>[...row.cells].map(cell=>Number(cell.textContent.replace(/−/g,'-')))));
      assert.deepEqual(target,certificate.target_B,id+' exact replay target');
    }
    await page.setViewportSize({width:390,height:844});
    await page.goto(root+'#r6-c88/exponents');await page.locator('#panel-exponents').scrollIntoViewIfNeeded();
    assert(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth));
    await page.screenshot({path:path.join(__dirname,'.qa/rank6-mobile-exponents.png')});
    assert.deepEqual(errors,[]);
    console.log('PASS: rank-six entries, three-digit navigation, proof, source links, 570-step Dynkin and infinite/surface replays, spectra and mobile layout.');
  } finally {await browser.close();}
})().catch(e=>{console.error(e);process.exitCode=1;});
