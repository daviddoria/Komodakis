// masks are defined as
// 0 = source region (NOT to be filled)
// 1 = target region (to be filled)


#ifndef Komodakis_H
#define Komodakis_H

#include "Types.h"

#include <itkImage.h>

struct Tile
{
  enum TileTypeEnum {BOUNDARY, SOURCE, TARGET};
  TileTypeEnum TileType;
  itk::ImageRegion<2> Region;
  unsigned int Id;
};

struct Node
{
  itk::Index<2> Location;
  unsigned int Id;
};

template <typename TImageType>
class Komodakis
{
public:
  Komodakis();

  void SetImage(typename TImageType::Pointer image);
  void SetMask(MaskType::Pointer mask);

  void Initialize();
  void CreateTiles();
  void CreateNodes();

  double SingleNodePotential(Tile tile, Node node); // SSD of single tile overlap with boundary
  double PairwisePotential(Tile tile1, Tile tile2, Node node1, Node node2); // SSD of overlapping region of two tiles

  bool IntersectTargetRegion(itk::ImageRegion<2> region);
  itk::ImageRegion<2> GetCenteredRegion(itk::Index<2> centerPixel, unsigned int edgeLength);

  std::vector<Node> GetNodes();
  std::vector<Tile> GetTiles();

private:
  typename TImageType::Pointer Image;
  typename MaskType::Pointer Mask;

  unsigned int TileEdgeLength;
  unsigned int Gap;

  std::vector<Tile> Tiles;
  std::vector<Node> Nodes;
};

#endif