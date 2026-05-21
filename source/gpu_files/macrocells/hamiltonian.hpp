class Hamiltonian
{
public:
    virtual void measure(std::size_t mstep) = 0;
    virtual void flushCorrelations(hostCorrelations& cpuCorrelations, std::size_t mstep) = 0;
    // virtual void copyToFortran() = 0;
    virtual ~Hamiltonian() = default;

protected:
    Hamiltonian() = default;
};


