"""Simulates an example network of 450 cell receiving two kinds of exernal input as defined in the configuration file"""
import sys
from bmtk.simulator import bionet
from ateam.sim import cell_functions
from ateam.sim.cell_reports import save_morph_single
from ateam.sim.setup import SimManager

def run(config_file):
    conf = bionet.Config.from_json(config_file, validate=True)
    conf.build_env()

    graph = bionet.BioNetwork.from_config(conf)
    sim = bionet.BioSimulator.from_config(conf, network=graph)
    save_morph_single(sim, 'morph.csv') # save morphology
   # sim.run()
    bionet.nrn.quit_execution()


if __name__ == '__main__':
    if __file__ != sys.argv[-1]:
        run(sys.argv[-1])
    else:
        run('config.json')
