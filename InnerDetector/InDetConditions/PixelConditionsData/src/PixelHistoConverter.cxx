/*
   Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
 */

#include "PixelConditionsData/PixelHistoConverter.h"

#include "StoreGate/StoreGateSvc.h"
#include "GaudiKernel/ISvcLocator.h"

#include "TAxis.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

#include <iostream>

PixelHistoConverter::PixelHistoConverter()
{
}

StatusCode PixelHistoConverter::SetHisto1D(const TH1* histo) {

  if (!histo) {
    return StatusCode::FAILURE;
  }

  if (!SetAxis(m_xAxis, histo->GetXaxis())) {
    return StatusCode::FAILURE;
  }

  /// fill the content
  const std::size_t xSize = m_xAxis.nBins;
  m_content.resize(xSize);
 
  for (std::size_t x = 0; x < xSize; ++x) {
    m_content.at(x) = histo->GetBinContent(x+1);
  }

  return StatusCode::SUCCESS;
}

StatusCode PixelHistoConverter::SetHisto2D(const TH2* histo) {
  if (!histo) {
    return StatusCode::FAILURE;
  }

  if (!SetAxis(m_xAxis, histo->GetXaxis())) {
    return StatusCode::FAILURE;
  }
  if (!SetAxis(m_yAxis, histo->GetYaxis())) {
    return StatusCode::FAILURE;
  }

  /// fill the content use linearized version for performance reasons
  const std::size_t xSize = m_xAxis.nBins;
  const std::size_t ySize = m_yAxis.nBins;
  m_content.resize(xSize*ySize);
 
  for (std::size_t x = 0; x < xSize; ++x) {
    for (std::size_t y = 0; y < ySize; ++y) {
      const std::size_t position = x + y*xSize;
      m_content.at(position) = histo->GetBinContent(x+1,y+1);
    }
  }

  return StatusCode::SUCCESS;
}

StatusCode PixelHistoConverter::SetHisto3D(const TH3* histo) {
  if (!histo) {
    return StatusCode::FAILURE;
  }

  if (!SetAxis(m_xAxis, histo->GetXaxis())) {
    return StatusCode::FAILURE;
  }
  if (!SetAxis(m_yAxis, histo->GetYaxis())) {
    return StatusCode::FAILURE;
  }
  if (!SetAxis(m_zAxis, histo->GetZaxis())) {
    return StatusCode::FAILURE;
  }

  /// fill the content use linearized version for performance reasons
  const std::size_t xSize = m_xAxis.nBins;
  const std::size_t ySize = m_yAxis.nBins;
  const std::size_t zSize = m_zAxis.nBins;
  m_content.resize(xSize*ySize*zSize);
 
  for (std::size_t x = 0; x < xSize; ++x) {
    for (std::size_t y = 0; y < ySize; ++y) {
      for (std::size_t z = 0; z < zSize; ++z) {
        const std::size_t position = x + xSize*(y + (ySize * z));
        m_content.at(position) = histo->GetBinContent(x+1,y+1,z+1);
      }
    }
  }

  return StatusCode::SUCCESS;
}

float PixelHistoConverter::GetContent(const std::size_t& x) const {
  return m_content[x];
}

float PixelHistoConverter::GetContent(const std::size_t& x, const std::size_t& y) const {
  const std::size_t position = x + y*(m_xAxis.nBins);
  return m_content[position];
}

float PixelHistoConverter::GetContent(const std::size_t& x, const std::size_t& y, const std::size_t& z) const {
  const std::size_t position = x + m_xAxis.nBins*(y + (m_yAxis.nBins * z));
  return m_content[position];
}

bool PixelHistoConverter::IsOverflowZ(const float value) const {
  return (value >= m_zAxis.max) ? true : false;
}

bool PixelHistoConverter::IsFirstZ(const float value) const {
  return (GetBinZ(value) == 0);
}

float PixelHistoConverter::GetBinX(const float value) const {
  return FindBin(m_xAxis, value);
}
  
float PixelHistoConverter::GetBinY(const float value) const {
  return FindBin(m_yAxis, value);
}

float PixelHistoConverter::GetBinZ(const float value) const {
  return FindBin(m_zAxis, value);
}

bool PixelHistoConverter::SetAxis(Axis& axis, const TAxis* rootAxis) {

  if (!rootAxis) {
    return false;
  }

  axis.nBins = rootAxis->GetNbins();
  axis.min   = rootAxis->GetXmin();
  axis.max   = rootAxis->GetXmax();

  if (axis.nBins < 1) {
    return false;
  }

  /// check if the histogram has equidistant bins
  const float width = rootAxis->GetBinWidth(1);
  for (std::size_t ibin = 2; ibin <= axis.nBins; ++ibin) {
    /// use a threshold for imperfect binning
    if (std::abs(rootAxis->GetBinWidth(ibin) - width) > 0.01*width) {
      return false;
    }
  }

  /// storing as 1/width to avoid (slow) division in retrieving
  axis.width = 1.*axis.nBins/(axis.max - axis.min);

  return true;
}

std::size_t PixelHistoConverter::FindBin(const Axis& axis, const float value) const {
  if (value <= axis.min) return 0;
  if (value >= axis.max) return (axis.nBins - 1);

  return ((value - axis.min) * axis.width);
}
