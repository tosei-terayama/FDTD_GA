#ifndef PML_H_
#define PML_H_

class pml{
private:
    int p_j1, p_j2, p_k1, p_k2;
    
public:
    void set_point_1(int, int);
    void set_point_2(int, int);
    void set_point(int v_y1, int v_y2, int v_z1, int v_z2);

    int j1(void){ return p_j1; }
    int j2(void){ return p_j2; }
    int k1(void){ return p_k1; }
    int k2(void){ return p_k2; }
};

#endif /* PML_H_ */
