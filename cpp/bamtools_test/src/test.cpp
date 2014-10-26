#include <string>
#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
using namespace BamTools;

int main(int argc, chr* argv[]) {
std::vector<std::string> inputFilenames;
std::string outputFilename;

 inputFilenames = argv[1];
 outputFilename = argv[2];

//provide input/output names somehow

 BamTools::BamMultiReader reader;
if ( !reader.Open(inputFilenames) ) {
  std::cerr << "Could not open input BAM files." << std::endl;
  return 0;
}

const SamHeader header = reader.GetHeader();
const RefVector references = reader.GetReferenceData();

BamWriter writer;
if ( !writer.Open(outputFilename, header, references) ) {
  std::cerr << "Could not open output BAM file" << std::endl;
  return 0;
}

BamAlignment al;
while ( reader.GetNextAlignmentCore(al) ) {
  if ( al.MapQuality >=90 ) writer.SaveAlignment(al);
}

reader.Close();
writer.Close();
}
