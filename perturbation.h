#ifndef PERTURBATION_H_
#define PERTURBATION_H_

class perturbation{
private:
    double Lati, Longi, Alti;
    int P_r0, P_th0, P_phi0;   // Perturbation center //
    double N0;  // Max enhacement //
    double Sigma_r;
    double Sigma_h;

public:
    void set_geo(double, double, double);
    void set_center(int, int, int);
    void set_alpha(double);
    void set_sigma(double, double);

    double lati(void){ return Lati; }
    double longi(void) { return Longi; }
    double alti(void) { return Alti; }
    double r0(void) { return P_r0; }
    double th0(void) { return P_th0; }
    double phi0(void) { return P_phi0; }
    double alpha(void) { return N0; }
    double sig_r(void) { return Sigma_r; }
    double sig_h(void) { return Sigma_h; }
};
#endif /* PERTURBATION_H_ */