#include <iostream>

void b(void);
void c(void);
void d(void);

int main(void){

    for(int i = 0; i < 1000; i++){
        std::cout << "\r" << i;
        b();
    }

    std::cout << std::endl;

    for(int i = 0; i < 10000; i++){
        std::cout << "\r" << i;
        c();
    }

    std::cout << std::endl;

    for(int i = 0; i < 50000; i++){
        std::cout << "\r" << i;
        d();
    }
    
    std::cout << std::endl << "process complete." << std::endl;

    return 0;
}