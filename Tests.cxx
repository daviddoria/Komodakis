#include "Types.h"
#include "Komodakis.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"

void WriteImage(ImageType::Pointer, std::string Filename);
void WriteImage(MaskType::Pointer, std::string Filename);

int main(int argc, char* argv[])
{
  typedef itk::ImageFileReader<ImageType> ImageReaderType;
  typedef itk::ImageFileReader<MaskType> MaskReaderType;

  std::string imageFilename = argv[1];
  std::cout << "Reading " << imageFilename << std::endl;

  ImageReaderType::Pointer imageReader = ImageReaderType::New();
  imageReader->SetFileName(imageFilename);
  imageReader->Update();

  WriteImage(imageReader->GetOutput(), "ImageRead.png");

  std::string maskFilename = argv[2];
  std::cout << "Reading " << maskFilename << std::endl;

  MaskReaderType::Pointer maskReader = MaskReaderType::New();
  maskReader->SetFileName(maskFilename);
  maskReader->Update();

  WriteImage(maskReader->GetOutput(), "MaskRead.png");

  Komodakis<ImageType> KomodakisMethod;
  KomodakisMethod.SetImage(imageReader->GetOutput());
  KomodakisMethod.SetMask(maskReader->GetOutput());
  KomodakisMethod.Initialize();

  std::vector<Node> nodes = KomodakisMethod.GetNodes();
  std::cout << "There are " << nodes.size() << " nodes." << std::endl;

  ImageType::Pointer nodeImage = ImageType::New();
  //nodeImage->Graft(imageReader->GetOutput());
  nodeImage->Graft(maskReader->GetOutput());

  for(unsigned int i = 0; i < nodes.size(); i++)
    {
    nodeImage->SetPixel(nodes[i].ImageLocation, 100);
    }

  WriteImage(nodeImage, "Nodes.png");

  KomodakisMethod.CompleteImage();

  std::cout << "PairwisePotentialCallCount: " << KomodakisMethod.PairwisePotentialCallCount << std::endl;
  std::cout << "SingleNodePotentialCallCount: " << KomodakisMethod.SingleNodePotentialCallCount << std::endl;

  for(unsigned int i = 0; i < nodes.size(); i++)
    {
    std::cout << "Node " << i << " : " << nodes[i].Label << std::endl;
    }

  return EXIT_SUCCESS;
}

void WriteImage(ImageType::Pointer image, std::string Filename)
{
  typedef  itk::ImageFileWriter< ImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(Filename);
  writer->SetInput(image);
  writer->Update();
}
/*
void WriteImage(MaskType::Pointer image, std::string Filename)
{
  typedef  itk::ImageFileWriter< MaskType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(Filename);
  writer->SetInput(image);
  writer->Update();
}
*/