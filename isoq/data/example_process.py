import isoq.process as iq

# data files
datafile = "example_dataset.csv"
calibfile = "example_calibration.csv"
folder = "D:/GIT/IsoQ/IsoQ/isoq/data/"

# processing parameters
tracer = '13C'
resolution = 70000
mz_of_resolution = 400
tracer_purity = None
correct_NA_tracer = False
resolution_formula_code = 'orbitrap'
purity15N = [0.01, 0.99]

# run processing
iq.run(folder, datafile, calibfile, tracer, resolution, mz_of_resolution, tracer_purity, correct_NA_tracer, resolution_formula_code, purity15N=purity15N, verbose=False)

