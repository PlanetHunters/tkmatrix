from tkmatrix.tkmatrix_class import MATRIX

ir = MATRIX("TIC 85400193", [18], "../sherlockpipe/run_tests/experiment/ir")\
    .plot_results(step_period=0.1, step_radius=0.1, phases=3, plot_step_period=5, plot_step_radius=5)
