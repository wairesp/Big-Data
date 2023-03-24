    #include <iostream>
    #include <fstream>
    #include <string>
    #include <unordered_map>
    #include <thread>
    #include <vector>
    #include <experimental/filesystem>
    #include <chrono>
    #include <random>

    void count_words(const std::string& chunk, std::unordered_map<std::string, int>& word_count) {
        std::istringstream iss(chunk);
        std::string word;
        while (iss >> word) {
            ++word_count[word];
        }
    }

    int main() {
        constexpr size_t chunk_size = 1024 * 1024; // 1 MB
        std::ifstream file("loremx2gb.txt");
        std::string file_path = "loremx2gb.txt";
        size_t file_size = std::experimental::filesystem::file_size(file_path);

        if (!file) {
            std::cerr << "Failed to open file\n";
            return -1;
        }

        std::vector<std::thread> threads;
        std::vector<std::unordered_map<std::string, int>> word_counts;
        word_counts.reserve(file_size / chunk_size + 1);

        std::random_device rd;  //Will be used to obtain a seed for the random number engine
        std::chrono::time_point<std::chrono::high_resolution_clock> tstart, end;
        tstart = std::chrono::high_resolution_clock::now();

        while (file) {
            auto chunk = std::string(chunk_size, '\0');
            file.read(&chunk[0], chunk_size);
            word_counts.emplace_back();
            threads.emplace_back(count_words, chunk, std::ref(word_counts.back()));
        }


        for (auto& thread : threads) {
            thread.join();
        }

        end = std::chrono::high_resolution_clock::now();
        int64_t duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - tstart).count();
        std::cout << "\nTiempo" << duration * 0.000000001f << " s." << "\n";

        std::ofstream outdata;
        outdata.open("output2GB.txt");
        outdata << "Tiempo: " << duration * 0.000000001f << " s." << "\n";

        //fill the map
        std::unordered_map<std::string, int> global_word_count;
        for (const auto& word_count : word_counts) {
            for (const auto& [word, count] : word_count) {
                global_word_count[word] += count;
            }
        }

        int x = global_word_count.size();
        std::cout << std::endl <<"Nro. Palabras:"<< x << std::endl;
        outdata << std::endl <<"Nro. Palabras: "<< x << std::endl << std::endl;

        for (auto const &pair: global_word_count) {
            std::cout << "{" << pair.first << ": " << pair.second << "}\n";
            outdata << "{" << pair.first << ": " << pair.second << "}\n";
        }

        outdata.close();
    }