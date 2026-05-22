class Hamiltonian
{
public:
    virtual void calculate(deviceLattice& gpuLattice, deviceEnergies& gpuEnergies, bool measure) = 0;
    // virtual void copyToFortran() = 0;
    virtual ~Hamiltonian() = default;

protected:
    Hamiltonian() = default;
};


