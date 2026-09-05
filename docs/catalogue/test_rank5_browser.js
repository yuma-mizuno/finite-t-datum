'use strict';
const assert=require('node:assert/strict'),path=require('node:path');
const {pathToFileURL}=require('node:url');
const {chromium}=require(process.env.PLAYWRIGHT_MODULE||'playwright');
(async()=>{
  const browser=await chromium.launch({headless:true,executablePath:process.env.BROWSER_EXECUTABLE});
  try {
    const page=await browser.newPage({viewport:{width:1450,height:1060}}),errors=[];
    page.on('pageerror',e=>errors.push(String(e)));
    const root=pathToFileURL(path.join(__dirname,'index.html')).href;
    await page.goto(root+'#r5-c45/matrices');
    assert.equal(await page.locator('.record-link').count(),55);
    await page.screenshot({path:path.join(__dirname,'.qa/rank5-matrices.png'),fullPage:true});
    await page.goto(root+'#r5-c54/exponents');
    await page.screenshot({path:path.join(__dirname,'.qa/rank5-exponents.png'),fullPage:true});
    await page.click('#tab-evidence');
    const links=await page.locator('#panel-evidence a').evaluateAll(xs=>xs.map(x=>x.getAttribute('href')));
    assert(!links.some(x=>x.includes('null')||x.includes('undefined')));
    assert(links.some(x=>x.endsWith('smt_queries.zip')));
    await page.click('#tab-proof');
    assert.match(await page.locator('#panel-proof').textContent(),/55 indecomposable families/);
    await page.screenshot({path:path.join(__dirname,'.qa/rank5-classification.png'),fullPage:true});
    await page.goto(root+'#r5-c54/notes');
    assert.match(await page.locator('#panel-notes').textContent(),/E8\^\(1,1\)/);
    await page.screenshot({path:path.join(__dirname,'.qa/rank5-notes.png'),fullPage:true});
    await page.setViewportSize({width:390,height:844});
    await page.goto(root+'#r5-c45/exponents');await page.locator('#panel-exponents').scrollIntoViewIfNeeded();
    assert(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth));
    await page.screenshot({path:path.join(__dirname,'.qa/rank5-mobile-exponents.png')});
    assert.deepEqual(errors,[]);console.log('PASS: rank-5 proof/source links, all family entries, notes, spectra and mobile layout.');
  } finally {await browser.close();}
})().catch(e=>{console.error(e);process.exitCode=1;});
