import pylhe
import numpy as np
import matplotlib.pyplot as plt

lhe_file = "/Users/albertodufour/MG5/mg5amcnlo/all_contributing/myLHE/unweighted_events.lhe"


def mll(p1, p2):
    p = np.array(p1) + np.array(p2)
    return np.sqrt(p[3]**2 - sum(p[l]**2 for l in range(3)))

def rap(p1, p2):
    p = np.array(p1) + np.array(p2)
    E  = p[3]
    pz = p[2]
    return 0.5 * np.log((E + pz) / (E - pz))

def cstar(p1, p2):
    # p1 = e-, p2 = e+
    p1 = np.array(p1)
    p2 = np.array(p2)
    
    # dilepton system
    p = p1 + p2
    E  = p[3]
    pz = p[2]
    mass = mll(p1, p2)  # or compute inline
    
    # boost the e- to the dilepton rest frame
    # boost vector along z
    beta  = pz / E
    gamma = E / mass
    
    # boosted e- momentum components
    pz1_boosted = gamma * (p1[2] - beta * p1[3])
    E1_boosted  = gamma * (p1[3] - beta * p1[2])
    
    # magnitude of 3-momentum of e- in rest frame
    p1_mag = np.sqrt(p1[0]**2 + p1[1]**2 + pz1_boosted**2)
    
    # direction of dilepton motion in lab = z axis
    # cos theta* is angle between e- and z in rest frame
    return pz1_boosted / p1_mag
    



if __name__ == "__main__":
    inv_mass = []
    rapidity = []
    cosstar = []
    for i, event in enumerate(pylhe.read_lhe_with_attributes(lhe_file)):
        # if i >= 100:
        #     break

        if not (i+1)%1000:
            print(f"\n--- Event {i+1} ---")
        final_state = [particle for particle in event.particles if particle.status == 1.0]
        pmu = [[final_state[k].px, final_state[k].py, final_state[k].pz, final_state[k].e] for k in range(2)]
        inv_mass.append(mll(pmu[0], pmu[1]))
        rapidity.append(rap(pmu[0], pmu[1]))
        cosstar.append(cstar(pmu[0], pmu[1]))
    
    plt.hist(inv_mass, bins=1000)
    plt.xlim(80, 100)
    plt.xlabel("m_ll")
    plt.savefig('/Users/albertodufour/code/DY2026/test/mll_dist_smeftsim_al.png')
    plt.clf()

    plt.hist(rapidity, bins=100)
    # plt.xlim(80, 100)l
    plt.xlabel("y")
    plt.savefig('/Users/albertodufour/code/DY2026/test/rap_dist_smeftNLO.png')
    plt.clf()

    plt.hist(cosstar, bins=100)
    # plt.xlim(80, 100)
    plt.xlabel("c*")
    plt.savefig('/Users/albertodufour/code/DY2026/test/cstar_dist_smeftNLO.png')

        