#ifndef GEOCOORDINATE_H_
#define GEOCOORDINATE_H_

class geocoordinate{
private:
    double Lati, Longi, Alti;
    int I, J, K;
    int Position_r, Position_th, Position_phi;

public:
    void geo_ijk(double, double, double);

    double lati(void){ return Lati; }
    double longi(void){ return Longi; }
    double alti(void){ return Alti; }
    int i(void){ return I; }
    int j(void){ return J; }
    int k(void){ return K; }
    int obs_r(void){ return Position_r; }
    int obs_th(void){ return Position_th; }
    int obs_phi(void){ return Position_phi; }

};

#endif /* GEOCOORDINATE_H_ */