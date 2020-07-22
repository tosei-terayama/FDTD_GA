#ifndef GEOCOORDINATE_H_
#define GEOCOORDINATE_H_

class geocoordinate{
private:
    double Lati, Longi, Alti;
    int Geo_I, Geo_J, Geo_K;
    int Obs_I, Obs_J, Obs_K;

public:
    void geo_ijk(double, double, double);
    void set_obs(int, int, int);

    double lati(void){ return Lati; }
    double longi(void){ return Longi; }
    double alti(void){ return Alti; }
    int geo_i(void){ return Geo_I; }
    int geo_j(void){ return Geo_J; }
    int geo_k(void){ return Geo_K; }
    int i(void){ return Obs_I; }
    int j(void){ return Obs_J; }
    int k(void){ return Obs_K; }

};

#endif /* GEOCOORDINATE_H_ */