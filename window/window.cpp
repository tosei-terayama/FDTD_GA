#include <windows.h>
#include <thread>
#include <chrono>

int WINAPI WinMain(HINSTANCE hInstance,     // 現在のコンストラクタ
                HINSTANCE hPrevInstance,       // 過去のコンストラクタ
                LPSTR lpCmdLine,            // コマンドライン
                int nCmdShow                // 表示状態
                ){

                    // ウインドウ作成 //
                    HWND hwnd = CreateWindow("STATIC", "WinAPITest", WS_CAPTION,
                                                0, 0,  //ウインドウの座標
                                                400, 300,  //ウインドウのサイズ
                                                NULL, NULL, hInstance, NULL);
                    if(hwnd == NULL) return 0;

                    // ウインドウ表示 //
                    ShowWindow(hwnd, SW_SHOW);

                    // 10秒スリープ //
                    std::this_thread::sleep_for(std::chrono::seconds(5));

                    return 0;
                }