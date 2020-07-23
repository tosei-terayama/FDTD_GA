#ifndef DATE_H_
#define DATE_H_

class date{
private:
    int Year, Month, Day;
    int MD;
    double Hour;

public:
    void set_y(int y);
    void set_m(int m);
    void set_d(int d);
    void set_h(double h);
    void set_ymd(int y, int m, int d);

    int iy(void){ return Year; }
    int im(void){ return Month; }
    int id(void){ return Day; }
    int imd(void){ return MD; }
    double ih(void){ return Hour; }

};

#endif /* DATE_H_ */