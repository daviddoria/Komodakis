#include "Komodakis.h"
#include "itkImage.h"

typedef itk::Image<unsigned char, 2> ImageType;

int main(int argc, char* argv[])
{
  Komodakis<ImageType> KomodakisMethod;

  return EXIT_SUCCESS;
}
