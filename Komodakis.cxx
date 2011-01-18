#include "Komodakis.h"

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkPasteImageFilter.h"
#include "itkOffset.h"
#include "itkImageFileWriter.h"

#include "mrf.h"
#include "MaxProdBP.h"

template <typename TImageType>
Komodakis<TImageType>::Komodakis()
{
  this->Image = TImageType::New();
  this->Mask = MaskType::New();

  this->TileEdgeLength = 10;
  this->Gap = this->TileEdgeLength/2;

  this->PairwisePotentialCallCount = 0;
  this->SingleNodePotentialCallCount = 0;
}

template <typename TImageType>
void Komodakis<TImageType>::SetImage(typename TImageType::Pointer image)
{
  this->Image->Graft(image);
}

template <typename TImageType>
void Komodakis<TImageType>::SetMask(MaskType::Pointer mask)
{
  this->Mask->Graft(mask);
}

template <typename TImageType>
std::vector<Node> Komodakis<TImageType>::GetNodes()
{
  return this->Nodes;
}

template <typename TImageType>
void Komodakis<TImageType>::CreateGridNodes()
{
  // find the upper left and lower right corner of nodes
  int minx = 1e5;
  int miny = 1e5;
  int maxx = -1e5;
  int maxy = -1e5;
  for(unsigned int i = 0; i < this->Nodes.size(); i++)
    {
    if(this->Nodes[i].ImageLocation[0] < minx)
      {
      minx = this->Nodes[i].ImageLocation[0];
      }
    if(this->Nodes[i].ImageLocation[1] < miny)
      {
      miny = this->Nodes[i].ImageLocation[1];
      }
    if(this->Nodes[i].ImageLocation[0] > maxx)
      {
      maxx = this->Nodes[i].ImageLocation[0];
      }
    if(this->Nodes[i].ImageLocation[1] > maxy)
      {
      maxy = this->Nodes[i].ImageLocation[1];
      }
    }
  std::cout << "Grid bounds: " << minx << " " << maxx << " " << miny << " " << maxy << std::endl;

  std::vector<Node> gridNodes;
  unsigned int nodeId = 0;
  unsigned int numXNodes = 0;
  for(int xid = 0; xid <= (maxx-minx)/this->Gap; xid++)
    {
    numXNodes++;
    for(int yid = 0; yid <= (maxy-miny)/this->Gap; yid++)
      {
      Node node;

      node.GridLocation[0] = xid;
      node.GridLocation[1] = yid;

      node.ImageLocation[0] = minx + (xid * this->Gap);
      node.ImageLocation[1] = miny + (yid * this->Gap);

      if(node.ImageLocation[0] < 4 || node.ImageLocation[1] < 4)
        {
        std::cout << "Strange" << std::endl;
        exit(-1);
        }
      node.Id = nodeId;

      gridNodes.push_back(node);
      nodeId++;
      }
    }

  this->GridSize[0] = numXNodes;
  this->GridSize[1] = gridNodes.size()/numXNodes;

  this->Nodes = gridNodes;
}

template <typename TImageType>
std::vector<Tile> Komodakis<TImageType>::GetTiles()
{
  return this->Tiles;
}

template <typename TImageType>
bool Komodakis<TImageType>::IntersectTargetRegion(itk::ImageRegion<2> region)
{
  itk::ImageRegionConstIterator<MaskType> imageIterator(this->Mask, region);

  unsigned int targetRegionCounter = 0;
  while(!imageIterator.IsAtEnd())
    {
    if(imageIterator.Get() != 0)
      {
      targetRegionCounter++;
      }

    ++imageIterator;
    }

  if(targetRegionCounter > 0)
    {
    return true;
    }
  else
    {
    return false;
    }
}

template <typename TImageType>
itk::ImageRegion<2> Komodakis<TImageType>::GetCenteredRegion(itk::Index<2> centerPixel, unsigned int edgeLength)
{
  itk::Size<2> size;
  size.Fill(this->TileEdgeLength);

  itk::Index<2> corner;
  corner[0] = centerPixel[0] - this->Gap;
  corner[1] = centerPixel[1] - this->Gap;

  return itk::ImageRegion<2>(corner, size);
}

template <typename TImageType>
void Komodakis<TImageType>::CreateTiles()
{
  unsigned int tileId = 0;

  //for(int i = this->Image->GetLargestPossibleRegion().GetIndex()[0]  + this->Gap; // don't start at the edge
  for(int i = this->Image->GetLargestPossibleRegion().GetIndex()[0]  + this->TileEdgeLength; // don't start at the edge
      i < this->Image->GetLargestPossibleRegion().GetIndex()[0] + static_cast<int>(this->Image->GetLargestPossibleRegion().GetSize()[0])
              //- static_cast<int>(this->Gap); // don't go all the way to the edge;
              - static_cast<int>(this->TileEdgeLength); // don't go all the way to the edge;
      i+=this->Gap)
    {
    //for(int j = this->Image->GetLargestPossibleRegion().GetIndex()[1]  + this->Gap; // don't start at the edge
    for(int j = this->Image->GetLargestPossibleRegion().GetIndex()[1]  + this->TileEdgeLength; // don't start at the edge
      j < this->Image->GetLargestPossibleRegion().GetIndex()[1] + static_cast<int>(this->Image->GetLargestPossibleRegion().GetSize()[1])
              //- static_cast<int>(this->Gap); // don't go all the way to the edge;
                - static_cast<int>(this->TileEdgeLength); // don't go all the way to the edge;
      j+=this->Gap)
      {

      itk::Size<2> tileSize;
      tileSize.Fill(this->TileEdgeLength);

      itk::Index<2> tileIndex; // the corner of the tile
      tileIndex[0] = i;
      tileIndex[1] = j;
      /*
      if(tileIndex[0] < 0 || tileIndex[1] < 0)
        {
        std::cerr << "Tile index is negative! This should not happen." << std::endl;
        exit(-1);
        }
      */
      itk::ImageRegion<2> tileRegion(tileIndex,tileSize);

      Tile tile;
      tile.Region = tileRegion;

      // Determine if tile is on the boundary
      itk::ImageRegionConstIterator<MaskType> maskIterator(this->Mask,tileRegion);

      unsigned int sourceRegionCounter = 0;
      unsigned int targetRegionCounter = 0;
      while(!maskIterator.IsAtEnd())
        {
        if(maskIterator.Get() == 0)
          {
          sourceRegionCounter++;
          }
        else
          {
          targetRegionCounter++;
          }

        ++maskIterator;
        }

      if(targetRegionCounter > 0 && sourceRegionCounter > 0)
        {
        tile.TileType = Tile::BOUNDARY;
        }
      else if(targetRegionCounter == 0)
        {
        tile.TileType = Tile::SOURCE;
        }
      else if(sourceRegionCounter == 0)
        {
        tile.TileType = Tile::TARGET;
        }
      else
        {
        std::cerr << "Invalid tile type!" << std::endl; // this should never happen
        exit(-1);
        }

      tile.Id = tileId;
      this->Tiles.push_back(tile);

      tileId++;
      }
    }

  std::cout << "There are " << this->Tiles.size() << " tiles." << std::endl;
}


template <typename TImageType>
void Komodakis<TImageType>::CreateNodes()
{
  unsigned int nodeId = 0;

  //for(int i = this->Image->GetLargestPossibleRegion().GetIndex()[0] + this->Gap; // don't start on the edge
  for(int i = this->Image->GetLargestPossibleRegion().GetIndex()[0] + this->TileEdgeLength; // don't start on the edge
      i < this->Image->GetLargestPossibleRegion().GetIndex()[0] + static_cast<int>(this->Image->GetLargestPossibleRegion().GetSize()[0])
                //- static_cast<int>(this->Gap); // don't go all the way to the edge
                - static_cast<int>(this->TileEdgeLength); // don't go all the way to the edge
      i+=this->Gap)
    {
    //for(int j = this->Image->GetLargestPossibleRegion().GetIndex()[1] + this->Gap; // don't start on the edge
    for(int j = this->Image->GetLargestPossibleRegion().GetIndex()[1] + this->TileEdgeLength; // don't start on the edge
      j < this->Image->GetLargestPossibleRegion().GetIndex()[1] + static_cast<int>(this->Image->GetLargestPossibleRegion().GetSize()[1])
                //- static_cast<int>(this->Gap); // don't go all the way to the edge
                  - static_cast<int>(this->TileEdgeLength); // don't go all the way to the edge
      j+=this->Gap)
      {

      itk::Index<2> nodeIndex;
      nodeIndex[0] = i;
      nodeIndex[1] = j;

      if(IntersectTargetRegion(GetCenteredRegion(nodeIndex, this->TileEdgeLength)))
        {
        Node node;

        node.Id = nodeId;
        node.ImageLocation = nodeIndex;
        this->Nodes.push_back(node);
        nodeId++;
        }

      }
    }
  std::cout << "There are " << this->Nodes.size() << " nodes." << std::endl;
}

template <typename TImageType>
double Komodakis<TImageType>::SingleNodePotential(Tile tile, Node node)
{
  //std::cout << "SingleNodePotential: Tile.Index = " << tile.Region.GetIndex() << std::endl;// << " Node = " << node << std::endl;
  //std::cout << "SingleNodePotential: Tile Region = " << tile.Region << std::endl;

  if(tile.TileType == TARGET)
    {
    return 0;
    }

  if(tile.TileType == SOURCE)
    {
    std::cerr << "Tile it outside of region - this should never happen
    exit(-1);
    }
  itk::ImageRegion<2> destination = GetCenteredRegion(node.ImageLocation, this->TileEdgeLength);

  if(destination.GetIndex()[0] < 0 || destination.GetIndex()[1] < 0)
    {
    std::cerr << "Destination index is negative! This should not happen." << std::endl;
    std::cerr << destination << std::endl;
    std::cerr << "Node id: " << node.Id << std::endl;
    std::cerr << "Node Image Location: " << node.ImageLocation << std::endl;
    std::cerr << "Node Grid Location: " << node.GridLocation << std::endl;
    exit(-1);
    }

  // SSD of single tile overlap with boundary
  itk::ImageRegionConstIterator<TImageType> tileImageIterator(this->Image, tile.Region);
  itk::ImageRegionConstIterator<TImageType> imageIterator(this->Image, destination);
  itk::ImageRegionConstIterator<MaskType> maskIterator(this->Mask, destination);

  double sum = 0;
  while(!imageIterator.IsAtEnd())
    {
    if(maskIterator.Get())
      {
      double difference = static_cast<double>(tileImageIterator.Get()) - static_cast<double>(imageIterator.Get());
      sum += difference * difference;
      }

    ++tileImageIterator;
    ++imageIterator;
    ++maskIterator;
    }

  this->SingleNodePotentialCallCount++;

  /*
  if(sum == 0)
    {
    std::cout << "Single: " << sum << " Node: " << node.Id << " Tile: " << tile.Id << std::endl;
    WriteTile(tile, node);
    //exit(-1);
    }
  */

  return sum;
}

template <typename TImageType>
double Komodakis<TImageType>::PairwisePotential(Tile tile1, Tile tile2, Node node1, Node node2)
{
  //std::cout << "PairwisePotential: Tile1.Index = " << tile1.Region.GetIndex() << " Tile2.Index = " << tile2.Region.GetIndex() << std::endl;
  // SSD of overlapping region of two tiles

  itk::ImageRegion<2> destination1 = GetCenteredRegion(node1.ImageLocation, this->TileEdgeLength);
  itk::ImageRegion<2> destination2 = GetCenteredRegion(node2.ImageLocation, this->TileEdgeLength);

  itk::Offset<2> offset1 = destination1.GetIndex() - tile1.Region.GetIndex();
  itk::Offset<2> offset2 = destination2.GetIndex() - tile2.Region.GetIndex();

  // Find the intersection of the two regions
  itk::ImageRegion<2> overlap = destination1;
  overlap.Crop(destination2);

  itk::ImageRegionConstIteratorWithIndex<TImageType> regionIterator(this->Image, overlap);
  itk::ImageRegionConstIteratorWithIndex<MaskType> maskIterator(this->Mask, overlap);

  double sum = 0;
  while(!regionIterator.IsAtEnd())
    {
    if(maskIterator.Get())
      {
      double difference = static_cast<double>(this->Image->GetPixel(regionIterator.GetIndex() + offset1)) -
                          static_cast<double>(this->Image->GetPixel(regionIterator.GetIndex() + offset2));
      sum += difference * difference;
      }

    ++regionIterator;
    ++maskIterator;
    }

  this->PairwisePotentialCallCount++;

  /*
  if(sum == 0)
    {
    std::cout << "Pairwise: " << sum << " Node1: " << node1.Id << " Node2: " << node2.Id << " Tile1: " << tile1.Id << " Tile2: " << tile2.Id << std::endl;
    WritePair(tile1, node1, tile2, node2);
    exit(-1);
    }
  */

  return sum;
}

template <typename TImageType>
void Komodakis<TImageType>::WritePair(Tile tile1, Node node1, Tile tile2, Node node2)
{
  // Create a black image
  ImageType::IndexType start;
  start.Fill(0);

  ImageType::RegionType region(start, this->Image->GetLargestPossibleRegion().GetSize());
  ImageType::Pointer blankImage = ImageType::New();
  blankImage->SetRegions(region);
  blankImage->Allocate();
  blankImage->FillBuffer(0);

  {
  typedef itk::PasteImageFilter <ImageType, ImageType >
    PasteImageFilterType;
  PasteImageFilterType::Pointer pasteFilter
    = PasteImageFilterType::New ();
  pasteFilter->SetSourceImage(this->Image);
  pasteFilter->SetDestinationImage(blankImage);
  pasteFilter->SetSourceRegion(tile1.Region);
  pasteFilter->SetDestinationIndex(node1.ImageLocation);
  pasteFilter->Update();

  typedef  itk::ImageFileWriter< ImageType  > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName("Tile1.png");
  writer->SetInput(pasteFilter->GetOutput());
  writer->Update();
  }

  {
  typedef itk::PasteImageFilter <ImageType, ImageType >
    PasteImageFilterType;
  PasteImageFilterType::Pointer pasteFilter
    = PasteImageFilterType::New ();
  pasteFilter->SetSourceImage(this->Image);
  pasteFilter->SetDestinationImage(blankImage);
  pasteFilter->SetSourceRegion(tile2.Region);
  pasteFilter->SetDestinationIndex(node2.ImageLocation);
  pasteFilter->Update();

  typedef  itk::ImageFileWriter< ImageType  > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName("Tile2.png");
  writer->SetInput(pasteFilter->GetOutput());
  writer->Update();
  }
}

template <typename TImageType>
void Komodakis<TImageType>::WriteTile(Tile tile, Node node)
{
  // Create a black image
  ImageType::IndexType start;
  start.Fill(0);

  ImageType::RegionType region(start, this->Image->GetLargestPossibleRegion().GetSize());
  ImageType::Pointer blankImage = ImageType::New();
  blankImage->SetRegions(region);
  blankImage->Allocate();
  blankImage->FillBuffer(0);

  {
  typedef itk::PasteImageFilter <ImageType, ImageType >
    PasteImageFilterType;
  PasteImageFilterType::Pointer pasteFilter
    = PasteImageFilterType::New ();
  pasteFilter->SetSourceImage(this->Image);
  pasteFilter->SetDestinationImage(blankImage);
  pasteFilter->SetSourceRegion(tile.Region);
  pasteFilter->SetDestinationIndex(node.ImageLocation);
  pasteFilter->Update();

  typedef  itk::ImageFileWriter< ImageType  > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName("TileOnBlank.png");
  writer->SetInput(pasteFilter->GetOutput());
  writer->Update();
  }

  {
  typedef itk::PasteImageFilter <ImageType, ImageType >
    PasteImageFilterType;
  PasteImageFilterType::Pointer pasteFilter
    = PasteImageFilterType::New ();
  pasteFilter->SetSourceImage(this->Image);
  pasteFilter->SetDestinationImage(this->Image);
  pasteFilter->SetSourceRegion(tile.Region);
  pasteFilter->SetDestinationIndex(node.ImageLocation);
  pasteFilter->Update();

  typedef  itk::ImageFileWriter< ImageType  > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName("TileOnImage.png");
  writer->SetInput(pasteFilter->GetOutput());
  writer->Update();
  }
}

template <typename TImageType>
void Komodakis<TImageType>::Initialize()
{
  CreateTiles();
  CreateNodes();
  CreateGridNodes();
}

template <typename TImageType>
MRF::CostVal Komodakis<TImageType>::dataCost(int node, int label)
{
  return SingleNodePotential(this->Tiles[label], this->Nodes[node]);
}

template <typename TImageType>
MRF::CostVal Komodakis<TImageType>::smoothnessCost(int node1, int node2, int label1, int label2)
{
  // The cost of assigning label1 to node1 and label2 to node2

  return PairwisePotential(this->Tiles[label1], this->Tiles[label2],
                           this->Nodes[node1], this->Nodes[node2]);
}

Komodakis<ImageType>* globalObject;

template <typename TImageType>
MRF::CostVal Komodakis<TImageType>::staticSmoothnessCost(int node1, int node2, int label1, int label2)
{
  return globalObject->smoothnessCost(node1, node2, label1, label2);
}

template <typename TImageType>
MRF::CostVal Komodakis<TImageType>::staticDataCost(int node, int label)
{
  return globalObject->dataCost(node, label);
}

template <typename TImageType>
void Komodakis<TImageType>::CompleteImage()
{
  globalObject = this;

  DataCost *data         = new DataCost(Komodakis<TImageType>::staticDataCost);
  SmoothnessCost *smoothness = new SmoothnessCost(Komodakis<TImageType>::staticSmoothnessCost);
  EnergyFunction *energy = new EnergyFunction(data,smoothness);

  MRF* mrf = new MaxProdBP(this->GridSize[0],this->GridSize[1],this->Tiles.size(),energy);

  mrf->initialize();
  mrf->clearAnswer();

  MRF::EnergyVal E = mrf->totalEnergy();
  std::cout << "Energy at the Start= " << (float)E << " (" << (float)mrf->smoothnessEnergy()
                                       << ", " << (float)mrf->dataEnergy() << ")" << std::endl;;

  float tot_t = 0;
  float t;
  for (int iter=0; iter < 10; iter++)
  {
    mrf->optimize(1, t);

    E = mrf->totalEnergy();
    tot_t = tot_t + t ;
    std::cout << "energy = " << (float)E << " (" << tot_t << " secs)" << std::endl;
  }

  for(unsigned int i = 0; i < this->Nodes.size(); i++)
    {
    this->Nodes[i].Label = mrf->getLabel(i);
    }

  delete mrf;
}


// Explicit instantiations
template class Komodakis<itk::Image<unsigned char, 2> >;
