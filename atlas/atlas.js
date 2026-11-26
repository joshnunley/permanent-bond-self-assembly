{
  const {kind, stock, rows} = JSON.parse(document.getElementById('atlas-data').textContent);
  const byId = new Map(rows.map(row => [row.id, row]));
  const gallery = document.getElementById('shapes');
  const detail = document.getElementById('detail');
  const status = document.getElementById('status');
  const search = document.getElementById('search');
  const previous = document.getElementById('previous');
  const next = document.getElementById('next');
  const range = document.getElementById('range');
  const namespace = 'http://www.w3.org/2000/svg';
  const pageSize = 48;
  let pageNumber = 0;
  let selection = null;

  function element(tag, text, className) {
    const node = document.createElement(tag);
    if (text) node.textContent = text;
    if (className) node.className = className;
    return node;
  }

  function vertices([x, y]) {
    if (kind === 'squares') return [[x, -y], [x + 1, -y], [x + 1, 1 - y], [x, 1 - y]];
    const i = Math.floor(x / 3);
    const j = Math.floor(y / 3);
    const up = x - 3 * i === 1;
    const triangle = up ? [[i, j], [i + 1, j], [i, j + 1]] : [[i + 1, j], [i, j + 1], [i + 1, j + 1]];
    return triangle.map(([a, b]) => [a + b / 2, -Math.sqrt(3) * b / 2]);
  }

  function drawing(cells) {
    const polygons = cells.map(vertices);
    const xs = polygons.flat().map(([x]) => x);
    const ys = polygons.flat().map(([, y]) => y);
    const x = Math.min(...xs);
    const y = Math.min(...ys);
    const pad = .12;
    const picture = document.createElementNS(namespace, 'svg');
    picture.setAttribute('class', 'shape');
    picture.setAttribute('aria-hidden', 'true');
    picture.setAttribute('viewBox', `${x - pad} ${y - pad} ${Math.max(...xs) - x + 2 * pad} ${Math.max(...ys) - y + 2 * pad}`);
    for (const polygon of polygons) {
      const cell = document.createElementNS(namespace, 'polygon');
      cell.setAttribute('points', polygon.map(point => point.join(',')).join(' '));
      picture.append(cell);
    }
    return picture;
  }

  function filteredRows() {
    const word = search.value.trim().toUpperCase();
    return rows.filter(row => {
      if (!row.id.includes(word)) return false;
      switch (status.value) {
        case 'yes': return row.factorizations.length > 0;
        case 'no': return row.factorizations.length === 0;
        case 'exact1': return row.target_state_bounds?.[1] === 1;
        case 'exact2': return row.target_state_bounds?.[0] === 2 && row.target_state_bounds[1] === 2;
        default: return true;
      }
    });
  }

  function stateCount(bounds) {
    const [lower, upper] = bounds;
    if (lower !== upper) return `${lower}–${upper} states (bounds)`;
    return `Minimum ${lower} ${lower === 1 ? 'state' : 'states'}`;
  }

  function shapeCard(row) {
    const link = element('a', '', 'shape-link');
    link.href = `#${row.id}`;
    if (selection === row) link.setAttribute('aria-current', 'true');
    const count = row.factorizations.length;
    link.append(drawing(row.shape), element('strong', row.id),
      element('span', count ? `${count} maximal ${count === 1 ? 'set' : 'sets'}` : 'Impossible', count ? 'yes' : 'no'));
    if (row.target_state_bounds) link.append(element('span', stateCount(row.target_state_bounds), 'muted'));
    link.addEventListener('click', () => {
      if (matchMedia('(max-width: 760px)').matches) requestAnimationFrame(() => detail.scrollIntoView());
    });
    return link;
  }

  function showGallery() {
    const candidates = filteredRows();
    const start = pageNumber * pageSize;
    const visible = candidates.slice(start, start + pageSize);
    gallery.replaceChildren(...visible.map(shapeCard));
    if (!visible.length) gallery.append(element('p', 'No matching shapes.'));
    previous.disabled = pageNumber === 0;
    next.disabled = start + pageSize >= candidates.length;
    range.textContent = visible.length ? `${start + 1}–${start + visible.length} of ${candidates.length}` : '0 shapes';
  }

  function partLink(row, index) {
    const part = row.parts[index];
    const id = row.part_targets[index];
    const hand = row.part_hands[index];
    const label = hand ? `${id} · ${hand}` : id;
    const item = element(part.length < stock ? 'a' : 'span');
    if (part.length < stock) item.href = `stock-${String(part.length).padStart(2, '0')}.html#${id}`;
    item.setAttribute('aria-label', `Part ${index + 1}, ${part.length} cells, shape ${label}`);
    item.append(drawing(part), element('span', label));
    return item;
  }

  function drawFactorization(row, group, index) {
    const container = element('details');
    container.open = index === 0;
    const bounds = row.factorization_state_bounds?.[index];
    const cost = bounds ? ` · ${stateCount(bounds)}` : '';
    container.append(element('summary', `Set ${index + 1} · ${group.length} shape types${cost}`));
    const sizes = new Map();
    for (const part of group) {
      const length = row.parts[part].length;
      if (!sizes.has(length)) sizes.set(length, []);
      sizes.get(length).push(part);
    }
    for (const [size, parts] of sizes) {
      const line = element('div', '', 'size-row');
      const items = element('div', '', 'part-list');
      items.append(...parts.map(part => partLink(row, part)));
      line.append(element('h4', `${size} ${size === 1 ? 'cell' : 'cells'}`), items);
      container.append(line);
    }
    return container;
  }

  function drawWitness(row, recipe, states) {
    const n = kind === 'triangles' ? 3 : 4;
    const label = port => states === 1 ? String(port) : `${port < n ? 'I' : 'B'}${port % n}`;
    const maximal = row.factorizations.findIndex(group =>
      recipe.parts.length === group.length && recipe.parts.every(part => group.includes(part)));
    const heading = maximal < 0 ? 'Smaller set' : `Maximal set ${maximal + 1}`;
    const card = element('details');
    card.open = true;
    card.append(element('summary', `${heading} · ${stateCount([states, states])}`));
    if (states === 2) card.append(element('p', 'I is the state of an isolated cell. A cell switches to B once it has a bond.'));
    const pairs = recipe.ports.map(([a, b]) => `${label(a)} ↔ ${label(b)}`).join(', ');
    card.append(element('p', pairs ? `Allow ${pairs}.` : 'The rule has no binding pairs.'));
    const parts = element('div', '', 'part-list');
    parts.append(...recipe.parts.map(part => partLink(row, part)));
    card.append(parts);
    return card;
  }

  function showDetail() {
    const row = selection;
    const count = row.factorizations.length;
    const header = element('div', '', 'target-head');
    const labels = element('div');
    labels.append(element('h2', row.id), element('strong', count ? 'Guaranteed' : 'Impossible', count ? 'yes' : 'no'));
    if (row.target_state_bounds) labels.append(element('span', stateCount(row.target_state_bounds), 'muted'));
    header.append(drawing(row.shape), labels);
    detail.replaceChildren(header);
    if (!count) {
      detail.append(element('p', 'No compatible factorization exists at this stock.'));
      return;
    }
    detail.append(element('h3', `${count} maximal ${count === 1 ? 'factorization' : 'factorizations'}`),
      element('p', `${row.mirror_groups} ${row.mirror_groups === 1 ? 'set' : 'sets'} up to reflection of the whole factorization.`, 'muted'));
    const list = element('div', '', 'factorizations');
    list.append(...row.factorizations.map((group, i) => drawFactorization(row, group, i)));
    detail.append(list);
    const witnesses = element('div', '', 'factorizations witnesses');
    for (const recipe of row.one_state_catalogues || []) witnesses.append(drawWitness(row, recipe, 1));
    for (const recipe of row.two_state_witnesses || []) witnesses.append(drawWitness(row, recipe, 2));
    if (witnesses.children.length) detail.append(element('h3', 'Local rules'), witnesses);
  }

  function chooseFromHash() {
    selection = byId.get(window.location.hash.slice(1)) || rows.find(row => row.factorizations.length) || rows[0];
    const index = filteredRows().indexOf(selection);
    if (index >= 0) pageNumber = Math.floor(index / pageSize);
    showGallery();
    showDetail();
  }

  for (const input of [search, status]) input.addEventListener('input', () => { pageNumber = 0; showGallery(); });
  previous.addEventListener('click', () => { pageNumber--; showGallery(); });
  next.addEventListener('click', () => { pageNumber++; showGallery(); });
  window.addEventListener('hashchange', chooseFromHash);
  chooseFromHash();
}
