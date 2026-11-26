"""Publish the exact atlas data as linked pages that also work from disk."""

from argparse import ArgumentParser
from html import escape
import json
from pathlib import Path

from permanent_atlas import COUNTS, free, mirror


def page(title, prefix, body, data=None):
    script = '' if data is None else (
        '<script type="application/json" id="atlas-data">' +
        json.dumps(data, separators=(',', ':')) + '</script>\n' +
        f'<script defer src="{prefix}atlas.js"></script>'
    )
    return f'''<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{escape(title)} · Permanent assembly atlas</title>
  <link rel="stylesheet" href="{prefix}atlas.css">
</head>
<body>
<main class="page">
{body}
</main>
{script}
</body>
</html>
'''


def totals(rows):
    guaranteed = sum(bool(row['factorizations']) for row in rows)
    return (len(rows), guaranteed, len(rows) - guaranteed,
            sum(len(row['factorizations']) for row in rows), sum(row['mirror_groups'] for row in rows))


def state_counts(rows):
    return tuple(sum(row.get('target_state_bounds') == [q, q] for row in rows) for q in (1, 2))


def stock_page(kind, stock, rows, limit):
    total, guaranteed, without_guarantee, factorizations, mirror_groups = totals(rows)
    stock_links = ' '.join(f'<a href="stock-{m:02}.html"' +
                           (' aria-current="page"' if m == stock else '') + f'>{m}</a>'
                           for m in range(1, limit + 1))
    state_summary = state_filters = ''
    if 'target_state_bounds' in rows[0]:
        one, two = state_counts(rows)
        state_summary = f'  <p>{one:,} targets with minimum 1 state · {two:,} certified with minimum 2 states</p>'
        state_filters = '<option value="exact1">Minimum 1 state</option><option value="exact2">Minimum 2 states certified</option>'
    body = f'''
<nav class="trail" aria-label="Breadcrumb"><a href="../index.html">Atlas</a> / <a href="index.html">{kind.title()}</a> / Stock {stock}</nav>
<header>
  <h1>{stock} {kind}</h1>
  <p>{total:,} shapes · {guaranteed:,} guaranteed · {without_guarantee:,} without a guaranteeing rule</p>
  <p>{factorizations:,} maximal factorization sets · {mirror_groups:,} up to reflection</p>
{state_summary}
</header>
<nav class="stocks" aria-label="Stock sizes">{stock_links}</nav>
<div class="workspace">
  <section class="gallery" aria-labelledby="shape-list-title">
    <h2 id="shape-list-title">Shapes</h2>
    <div class="filters">
      <label>Result<select id="status"><option value="all">All shapes</option><option value="yes">Guaranteed</option><option value="no">No guaranteeing rule</option>{state_filters}</select></label>
      <label>Shape ID<input id="search" type="search" autocomplete="off" placeholder="e.g. {rows[0]['id']}"></label>
    </div>
    <div class="pagination"><button type="button" id="previous">Previous</button><span id="range" aria-live="polite"></span><button type="button" id="next">Next</button></div>
    <div id="shapes" class="tiles"></div>
  </section>
  <section id="detail" class="detail" aria-label="Selected shape"></section>
</div>
<noscript><p>Enable JavaScript to browse shapes. The exact shape data is also embedded in this page.</p></noscript>
'''
    title = f'Stock {stock} · {kind.title()}'
    return page(title, '../', body, {'kind': kind, 'stock': stock, 'rows': rows})


def publish(data, output):
    links = []
    has_states = False
    for kind in COUNTS:
        folders = sorted((data / kind).glob('[0-9][0-9].json'))
        if not folders:
            continue
        stocks = {int(path.stem): json.loads(path.read_text()) for path in folders}
        with_states = all('target_state_bounds' in rows[0] for rows in stocks.values())
        has_states |= with_states
        lookup = {tuple(map(tuple, row['shape'])): row['id'] for rows in stocks.values() for row in rows}
        destination = output / kind
        destination.mkdir(parents=True, exist_ok=True)
        for stock, rows in stocks.items():
            for row in rows:
                parts = [tuple(map(tuple, part)) for part in row['parts']]
                row['part_targets'] = [lookup[free(kind, part)] for part in parts]
                row['part_hands'] = ['' if part == mirror(kind, part) else ('A' if part < mirror(kind, part) else 'B')
                                     for part in parts]
            (destination / f'stock-{stock:02}.html').write_text(stock_page(kind, stock, rows, max(stocks)))
        table = []
        for stock, rows in stocks.items():
            counts = totals(rows) + (state_counts(rows) if with_states else ())
            figures = ''.join(f'<td>{n:,}</td>' for n in counts)
            table.append(f'<tr><th scope="row"><a href="stock-{stock:02}.html">{stock}</a></th>{figures}</tr>')
        state_columns = '<th scope="col">Minimum 1 state</th><th scope="col">Minimum 2 states certified</th>' if with_states else ''
        title = f'{kind.title()} at stocks 1–{max(stocks)}'
        body = f'''
<nav class="trail" aria-label="Breadcrumb"><a href="../index.html">Atlas</a> / {kind.title()}</nav>
<header><h1>{title}</h1><p>Guaranteed targets and all inclusion-maximal sets of permanent intermediate shapes.</p></header>
<div class="table-wrap"><table>
<thead><tr><th scope="col">Stock</th><th scope="col">Shapes</th><th scope="col">Guaranteed</th><th scope="col">No guaranteeing rule</th><th scope="col">Maximal sets</th><th scope="col">Up to reflection</th>{state_columns}</tr></thead>
<tbody>{''.join(table)}</tbody>
</table></div>
'''
        (destination / 'index.html').write_text(page(title, '../', body))
        links.append(f'<li><a href="{kind}/index.html">{title}</a></li>')
    state_intro = ''
    if has_states:
        state_intro = '''
<p>Local state counts use numbered edges and let neighboring cells settle before the next join. Each guaranteed target and maximal set has an exact minimum or an interval containing it. The target minimum allows any successful set of intermediates. Compare the two rules for the <a href="squares/stock-04.html#S04-0005">four-square target</a>.</p>
'''
    body = f'''
<header><h1>Permanent assembly atlas</h1><p>Which shapes are guaranteed from exactly their own stock of identical loose cells?</p></header>
<nav aria-label="Polygon atlases"><ul>{''.join(links)}</ul></nav>
<h2>How to read a shape page</h2>
<p>Each stock page lists all edge-connected target shapes up to rotation, translation, and reflection. Open a guaranteed target to see all its maximal sets of permitted parts. Intermediate parts keep their proper orientation class, so a reflected part can occupy a separate entry. The factorization count up to reflection identifies a set with its whole mirror image.</p>
<p>A shape listed under “No guaranteeing rule” can still be produced in some histories. <a href="../proofs/fixed-target.md#possibility-guarantee-and-probability">How guarantees differ from success probabilities</a></p>
<p>A smaller part links to the page where it is treated as a target in its own right. A part may occur inside a successful large-stock rule even when it cannot itself be guaranteed at its own stock.</p>
{state_intro}
<p><a href="README.md">Method, verification, and reproduction</a></p>
'''
    (output / 'index.html').write_text(page('Permanent assembly atlas', '', body))


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('data', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    publish(args.data, args.output)
