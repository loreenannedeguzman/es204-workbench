"""Plot presentation only: never changes the numerical tables."""
import math


def draw_chart(axis, chart, error_mode='Signed', xscale='Auto', yscale='Auto', tolerance=True):
    is_error = chart.get('error_plot', False)
    signed = error_mode == 'Signed'
    xs = chart.get('xscale', 'linear') if xscale == 'Auto' else xscale.lower()
    ys = chart.get('yscale', 'linear') if yscale == 'Auto' else yscale.lower()
    if is_error and yscale == 'Auto':
        ys = 'symlog' if signed else 'log'
    curves = chart['curves']
    values = [abs(y) for c in curves for y in c['y'] if y and math.isfinite(y)]
    threshold = min(values) if values else 1e-6
    tol = chart.get('tolerance')
    if tol and tol > 0:
        threshold = min(threshold, tol / 10)
    axis.set_xscale(xs)
    axis.set_yscale(ys, **({'linthresh': threshold} if ys == 'symlog' else {}))
    markers = ['s', 'x', 'o', '^', 'D', '+']
    lines = ['-', '--', '-.', ':']
    omitted = 0
    for i, c in enumerate(curves):
        pairs = []
        for x, y in zip(c['x'], c['y']):
            if is_error and not signed:
                y = abs(y)
            if (xs == 'log' and x <= 0) or (ys == 'log' and y <= 0):
                omitted += 1
                pairs.append((x, float('nan')))
            else:
                pairs.append((x, y))
        style = c.get('style', 'line')
        axis.plot([p[0] for p in pairs], [p[1] for p in pairs], label=c['label'],
                  linestyle='None' if style == 'points' else lines[i % len(lines)],
                  marker=markers[i % len(markers)] if style in ('nodes', 'points') else None,
                  markersize=5)
    if is_error and tolerance and tol and tol > 0:
        axis.axhline(tol, color='gray', linestyle=':', label=f'Tolerance {tol:g}%')
        if signed and ys != 'log':
            axis.axhline(-tol, color='gray', linestyle=':')
    ylabel = ('Signed' if signed else 'Absolute') + ' true error (%)' if is_error else chart['ylabel']
    axis.set(title=chart['title'], xlabel=chart.get('xlabel', 'x'), ylabel=ylabel)
    axis.grid(True, which='both', alpha=.3)
    axis.legend(fontsize=8)
    if omitted:
        axis.text(.01, .01, 'Log scale omits nonpositive values; use Symlog to show them.',
                  transform=axis.transAxes, fontsize=8, bbox=dict(facecolor='white', alpha=.9))
    return axis
