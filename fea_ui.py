"""Inputs for a general, two-node, one-dimensional Galerkin workflow."""
from tkinter import ttk


def build(app, ScrollForm):
    form = ScrollForm(app.tabs)
    app.tabs.add(form, text='FEA')
    p = form.body
    app.hint(p, "1D scalar, linear, two-node Galerkin elements. First assemble and inspect the working. Then manually partition the matrix and enter the reduced system below. No assignment-specific values are supplied.")
    app.field(p, 'fea_form', 'Equation form', 'Self-adjoint', choices=['Self-adjoint', 'Standard'])
    app.hint(p, "Self-adjoint: -(p*y\u2032)\u2032 + r*y = f. Standard: a*y\u2033 + b*y\u2032 + c*y = f; the workbench derives an integrating factor. Enter expressions in x; ln(x) or log(x), exp(x), sin(x), cos(x) and ^ are supported.")
    for key, label, value in [('p','Self-adjoint p(x)','1'), ('r','Self-adjoint r(x)','0'), ('a','Standard a(x)','1'), ('b','Standard b(x)','0'), ('c','Standard c(x)','0'), ('f','Right side f(x)','0')]:
        app.field(p, 'fea_'+key, label, value)
    app.field(p, 'fea_nodes', 'Increasing node coordinates, separated by commas', '')
    app.text(p, 'fea_regions', 'Optional self-adjoint regions: left; right; p(x); r(x); f(x), one per line. Must cover the domain; put every interface in the mesh.', 4)
    app.text(p, 'fea_loads', 'Optional forcing regions (either form): left; right; f(x). Must cover the domain. Cannot combine with the region table above.', 4)
    app.field(p, 'fea_integration', 'Element integration', 'Auto', choices=['Auto', 'Gauss 2', 'Analytical only'])
    app.hint(p, 'Auto uses supported symbolic antiderivatives and otherwise two-point Gauss quadrature. The report identifies the method for each integral. Nodes are numbered from 1; adjacent nodes form an element.')
    app.hint(p, "Endpoint conditions: alpha*y + beta*y\u2032 = gamma. The derivative is in increasing x at BOTH ends. For prescribed y, use beta=0. For prescribed y\u2032, use alpha=0.")
    for side in ('left','right'):
        for key, value in [('alpha','1'),('beta','0'),('gamma','0')]:
            app.field(p, 'fea_'+side+'_'+key, side.title()+' '+key, value)
    app.text(p, 'fea_point_loads', 'Optional additional weak-form nodal loads: node ID, load. Positive loads add to the assembled right side.', 3)
    app.field(p, 'fea_workflow', 'Workflow', 'Assemble only', choices=['Assemble only', 'Check partition and solve'])
    app.hint(p, 'After assembly: identify prescribed and free nodes. Form K_FF and F_F - K_FD*y_D yourself. Enter them below in your chosen free-node order. The workbench checks your partition before solving your entered system.')
    app.field(p, 'fea_free_order', 'Free node IDs, in your reduced-system order', '')
    app.text(p, 'fea_manual_matrix', 'Your reduced matrix: one row per line, spaces or commas', 5)
    app.text(p, 'fea_manual_rhs', 'Your reduced right side: spaces, commas or newlines', 3)
    app.field(p, 'fea_partition_tolerance', 'Absolute tolerance for checking entered matrix/RHS entries', '1e-6')
    ttk.Button(p, text='Run FEA workflow', command=app.run).pack(fill='x', pady=5)
    app.hint(p, 'Outputs: mesh, shapes, transformed equation, integration details, element matrices, assembled system and boundary terms. Solving adds nodal y, element slopes, residuals and plots. Use Rendered equations and Export equations for the working. Scope: regular 1D linear scalar equations; no nonlinear, 2D or higher-order elements.')
