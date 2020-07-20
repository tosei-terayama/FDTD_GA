#ifndef PERTURBATION_H_
#define PERTURBATION_H_

class perturbation{
private:
    double Lati, Longi, Alti;
    double P_r0, P_th0, P_phi0;   // Perturbation center //
    double N0;  // Max enhacement //
    double Sigma_r;
    double Sigma_h;
    int Range_r, Range_th, Range_phi;  // Perturbation range //

public:
    void set_geo(double, double, double);
    void set_center(double, double, double);
    void set_alpha(double);
    void set_sigma(double, double);
    void set_range(int, int, int);

    double lati(void){ return Lati; }
    double longi(void) { return Longi; }
    double alti(void) { return Alti; }
    double r0(void) { return P_r0; }
    double th0(void) { return P_th0; }
    double phi0(void) { return P_phi0; }
    double alpha(void) { return N0; }
    double sig_r(void) { return Sigma_r; }
    double sig_h(void) { return Sigma_h; }
    int range_r(void) { return Range_r; }
    int range_th(void) { return Range_th; }
    int range_phi(void) { return Range_phi; }
};
#endif /* PERTURBATION_H_ */