// Render report pages to PNG/PDF with headless Chromium.
//   NODE_PATH=$(npm root -g) node report/render.js shot <dist html> <out.png> <width> <height> [fullPage] [query]
//   NODE_PATH=$(npm root -g) node report/render.js screens   -> report/dist/screens/screen-1..5.png
//   NODE_PATH=$(npm root -g) node report/render.js poster    -> report/dist/poster.png and poster.pdf
const { chromium } = require("playwright");
const fs = require("fs"), path = require("path");
const root = path.resolve(__dirname);

function wrap(file) {
  // artifacts are wrapped in a document skeleton at publish time; do the same locally
  file = path.resolve(file); const body = fs.readFileSync(file, "utf8");
  const out = file.replace(/\.html$/, ".wrapped.html");
  fs.writeFileSync(out, `<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1,viewport-fit=cover"></head><body>${body}</body></html>`);
  return out;
}

async function launch() {
  const exe = fs.existsSync("/opt/pw-browsers/chromium") ? undefined : undefined;
  return chromium.launch({ args: ["--disable-web-security"], executablePath: exe });
}

async function shot(file, out, w, h, full, query) {
  const browser = await launch();
  const page = await browser.newPage({ viewport: { width: +w, height: +h }, deviceScaleFactor: 1 });
  const errs = [];
  page.on("console", m => { if (m.type() === "error") errs.push(m.text()); });
  page.on("pageerror", e => errs.push(String(e)));
  await page.goto("file://" + wrap(file) + (query || ""), { waitUntil: "networkidle" });
  await page.waitForTimeout(1500);
  await page.screenshot({ path: out, fullPage: full === "1" });
  const overflow = await page.evaluate(() => document.documentElement.scrollWidth - window.innerWidth);
  console.log(JSON.stringify({ out, overflow_px: overflow, errors: errs.slice(0, 5) }));
  await browser.close();
}

async function artboards(file, outdir, selector, scale, pdf) {
  const browser = await launch();
  const page = await browser.newPage({ viewport: { width: 2000, height: 2000 }, deviceScaleFactor: scale });
  const errs = [];
  page.on("pageerror", e => errs.push(String(e)));
  await page.goto("file://" + wrap(file) + "?still", { waitUntil: "networkidle" });
  await page.waitForTimeout(2000);
  fs.mkdirSync(outdir, { recursive: true });
  const boards = await page.$$(selector);
  for (let i = 0; i < boards.length; i++) {
    const id = await boards[i].getAttribute("data-name");
    const p = path.join(outdir, `${id || "board-" + (i + 1)}.png`);
    await boards[i].screenshot({ path: p });
    console.log("wrote", p);
  }
  if (pdf) {
    const box = await boards[0].boundingBox();
    await page.pdf({ path: pdf, width: `${box.width}px`, height: `${box.height}px`, printBackground: true, pageRanges: "1" });
    console.log("wrote", pdf);
  }
  if (errs.length) console.log("errors", errs);
  await browser.close();
}

(async () => {
  const [cmd, ...a] = process.argv.slice(2);
  if (cmd === "shot") await shot(...a);
  else if (cmd === "screens") await artboards(path.join(root, "dist", "screens.html"), path.join(root, "dist", "screens"), ".board", 1, null);
  else if (cmd === "poster") await artboards(path.join(root, "dist", "poster.html"), path.join(root, "dist"), ".poster", 2, path.join(root, "dist", "poster.pdf"));
})();
