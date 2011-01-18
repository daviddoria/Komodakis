// masks are defined as
// 0 = source region (NOT to be filled)
// 1 = target region (to be filled)

// This code is based on
// Image completion using global optimization
// Komodakis, N.

#ifndef Komodakis_H
#define Komodakis_H

#include "Types.h"

#include "mrf.h"

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
  itk::Index<2> GridLocation;
  itk::Index<2> ImageLocation;
  bool Boundary; // Is the node on the boundary of the region to complete?
  unsigned int Id;
  unsigned int Label;
};

template <typename TImageType>
class Komodakis
{
public:
  Komodakis();

  void SetImage(typename TImageType::Pointer image);
  void SetMask(MaskType::Pointer mask);

  void Initialize();

  std::vector<Node> GetNodes();
  std::vector<Tile> GetTiles();

  void CreateGridNodes();

  void CompleteImage();

  MRF::CostVal smoothnessCost(int node1, int node2, int label1, int label2);
  MRF::CostVal dataCost(int node, int label);

  static MRF::CostVal staticSmoothnessCost(int node1, int node2, int label1, int label2);
  static MRF::CostVal staticDataCost(int node, int label);


  unsigned int PairwisePotentialCallCount;
  unsigned int SingleNodePotentialCallCount;

private:
  typename TImageType::Pointer Image;
  typename MaskType::Pointer Mask;

  void WriteTile(Tile tile, Node node);
  void WritePair(Tile tile1, Node node1, Tile tile2, Node node2);

  void CreateTiles();
  void CreateNodes();

  double SingleNodePotential(Tile tile, Node node); // SSD of single tile overlap with boundary
  double PairwisePotential(Tile tile1, Tile tile2, Node node1, Node node2); // SSD of overlapping region of two tiles

  bool IntersectTargetRegion(itk::ImageRegion<2> region);
  itk::ImageRegion<2> GetCenteredRegion(itk::Index<2> centerPixel, unsigned int edgeLength);


  unsigned int TileEdgeLength;
  unsigned int Gap;

  std::vector<Tile> Tiles;
  std::vector<Node> Nodes;

  itk::Size<2> GridSize;

};

#endif