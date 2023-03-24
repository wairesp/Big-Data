#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <chrono>
#include <random>

int main() {
    std::unordered_map<std::string, int> wordCount;
    std::string word;
    std::ifstream file("loremx2gb.txt");

    std::ofstream outdata;
    outdata.open("outputSecuecial2GB.txt");


    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::chrono::time_point<std::chrono::high_resolution_clock> tstart, end;
    tstart = std::chrono::high_resolution_clock::now();


    while (file >> word) {
        ++wordCount[word];
    }

    end = std::chrono::high_resolution_clock::now();
    int64_t duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - tstart).count();
    std::cout << "Tiempo: " << duration * 0.000000001f << " s." << "\n";
    outdata << "Tiempo: " << duration * 0.000000001f << " s." << "\n";

    outdata << std::endl <<"Nro. Palabras: "<< wordCount.size() << std::endl << std::endl;

    for (const auto &pair : wordCount) {
        std::cout << "{" << pair.first << ": " << pair.second << "}\n";
        outdata << "{" << pair.first << ": " << pair.second << "}\n";
    }



    return 0;
}