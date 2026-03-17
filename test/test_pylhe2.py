import pylhe
import numpy as np

lhe_file = "/Users/albertodufour/MG5/mg5amcnlo/DrellYan/myLHEdir/DY_test_unweighted_events.lhe"


def mll(p1, p2):
    p = np.array(p1) + np.array(p2)
    return np.sqrt(p[3]**2 - sum(p[l]**2 for l in range(3)))

def rapidity(p1, p2):
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
    for i, event in enumerate(pylhe.read_lhe_with_attributes(lhe_file)):
        if i >= 3:
            break

       
        print(f"\n--- Event {i+1} ---")
        final_state = [particle for particle in event.particles if particle.status == 1.0]
        pmu = [[final_state[k].px, final_state[k].py, final_state[k].pz, final_state[k].e] for k in range(2)]

        print(f"mll= {mll(pmu[0], pmu[1])}")
        print(f"rapidity= {rapidity(pmu[0], pmu[1])}")
        print(f"cstar= {cstar(pmu[0], pmu[1])}")


        