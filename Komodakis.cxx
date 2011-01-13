#include "Komodakis.h"

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"

#include "itkOffset.h"

template <typename TImageType>
Komodakis<TImageType>::Komodakis()
{
  this->Image = TImageType::New();
  this->Mask = MaskType::New();

  this->TileEdgeLength = 10;
  this->Gap = this->TileEdgeLength/2;
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

  for(int i = this->Image->GetLargestPossibleRegion().GetIndex()[0];
      i < this->Image->GetLargestPossibleRegion().GetIndex()[0] + static_cast<int>(this->Image->GetLargestPossibleRegion().GetSize()[0])
              - static_cast<int>(this->Gap); // don't go all the way to the edge;
      i+=this->Gap)
    {
    for(int j = this->Image->GetLargestPossibleRegion().GetIndex()[1];
      j < this->Image->GetLargestPossibleRegion().GetIndex()[1] + static_cast<int>(this->Image->GetLargestPossibleRegion().GetSize()[1])
              - static_cast<int>(this->Gap); // don't go all the way to the edge;
      j+=this->Gap)
      {

      itk::Size<2> tileSize;
      tileSize.Fill(this->TileEdgeLength);

      itk::Index<2> tileIndex; // the corner of the tile
      tileIndex[0] = i;
      tileIndex[1] = j;

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
}


template <typename TImageType>
void Komodakis<TImageType>::CreateNodes()
{
  unsigned int nodeId = 0;

  for(int i = this->Image->GetLargestPossibleRegion().GetIndex()[0] + this->Gap; // don't start on the edge
      i < this->Image->GetLargestPossibleRegion().GetIndex()[0] + static_cast<int>(this->Image->GetLargestPossibleRegion().GetSize()[0])
                - static_cast<int>(this->Gap); // don't go all the way to the edge
      i+=this->Gap)
    {
    for(int j = this->Image->GetLargestPossibleRegion().GetIndex()[1] + this->Gap; // don't start on the edge
      j < this->Image->GetLargestPossibleRegion().GetIndex()[1] + static_cast<int>(this->Image->GetLargestPossibleRegion().GetSize()[1])
                - static_cast<int>(this->Gap); // don't go all the way to the edge
      j+=this->Gap)
      {

      itk::Index<2> nodeIndex;
      nodeIndex[0] = i;
      nodeIndex[1] = j;

      if(IntersectTargetRegion(GetCenteredRegion(nodeIndex, this->TileEdgeLength)))
        {
        Node node;

        node.Id = nodeId;
        node.Location = nodeIndex;
        this->Nodes.push_back(node);
        nodeId++;
        }

      }
    }
}

template <typename TImageType>
double Komodakis<TImageType>::SingleNodePotential(Tile tile, Node node)
{
  itk::ImageRegion<2> destination = GetCenteredRegion(node.Location, this->TileEdgeLength);

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

  return sum;
}

template <typename TImageType>
double Komodakis<TImageType>::PairwisePotential(Tile tile1, Tile tile2, Node node1, Node node2)
{
  // SSD of overlapping region of two tiles

  itk::ImageRegion<2> destination1 = GetCenteredRegion(node1.Location, this->TileEdgeLength);
  itk::ImageRegion<2> destination2 = GetCenteredRegion(node2.Location, this->TileEdgeLength);

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

  return sum;
}

template <typename TImageType>
void Komodakis<TImageType>::Initialize()
{
  CreateTiles();
  CreateNodes();
}

// Explicit instantiations
template class Komodakis<itk::Image<unsigned char, 2> >;
