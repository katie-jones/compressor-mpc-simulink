#include <sstream>
#include <fstream>
#include <iostream>
#include <string>

void scale_plot(std::ifstream& infile, std::ofstream& outfile,
                const double tmin, const double tmax) {
  std::string line;
  double t;
  const std::string start_plot = "addplot";
  std::string temp = "";
  bool points_found = false;

  // find start of plotting
  while (std::getline(infile, line)) {
    if (line.find(start_plot) != std::string::npos) {
      break;
    }
    outfile << line << std::endl;
  }

  // print any non-number lines after start of plotting
  do {
    std::istringstream iss(line);
    if (!(iss >> t)) {
      temp += line + "\n";
    } else {
      // if (t <= tmax && t >= tmin) {
        // points_found = true;
        // outfile << temp;
        // outfile << line << std::endl;
      // }
      break;
    }
  } while (std::getline(infile, line));

  // Read in times until the end
  do {
    std::istringstream iss(line);
    if (!(iss >> t)) {
      if (points_found) {
        outfile << line << std::endl;
      }
      break;
    }  // error
    if (t <= tmax && t >= tmin) {
      if (!points_found) {
        points_found = true;
        outfile << temp;
      }
      outfile << line << std::endl;
    }
  } while (std::getline(infile, line));
}

int main(int argc, char** argv) {
  if (argc < 4) {
    std::cerr << "Filename and time bounds must be specified." << std::endl;
    return 1;
  }

  const std::string filename(argv[1]);
  double tmin, tmax;

  std::istringstream ss(argv[2]);
  if (!(ss >> tmin)) std::cerr << "Invalid number " << argv[1] << '\n';
  std::istringstream ss2(argv[3]);
  if (!(ss2 >> tmax)) std::cerr << "Invalid number " << argv[1] << '\n';

  std::cout << "Bounds are: " << tmin << ", " << tmax << std::endl;

  std::ifstream infile;
  std::ofstream outfile;
  std::string output_filename = filename + ".new";
  infile.open(filename.c_str());
  outfile.open(output_filename.c_str());

  // print rest of file
  while (!(infile.eof() || infile.fail())) {
    scale_plot(infile, outfile, tmin, tmax);
  }

  return 0;
}
