#include "Config.H"

#include <exception>
#include <iostream>

int main(int argc, char** argv)
{
    if (argc != 2) {
        std::cerr << "usage: mhd2d_config_validate <config.json>\n";
        return 2;
    }
    try {
        const auto config = mhd::SimConfig::from_json_file(argv[1]);
        std::cout << "config valid: " << config.problem << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
