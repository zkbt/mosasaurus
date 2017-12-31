# how should spectra be shifted relative to each other?
correlationAnchors 8498.0 8542.0 8662.0
correlationRange   8350 8800
correlationSmooth  2

self.correlationAnchors = [float(x) for x in dictionary['correlationAnchors']]
self.correlationRange = [float(x) for x in dictionary['correlationRange']]
self.correlationSmooth = float(dictionary['correlationSmooth'])




# COSMIC RAY MITIGATION
self.cosmicThreshold = float(dictionary['cosmicThreshold'])
self.cosmicAbandon = float(dictionary['cosmicAbandon'])
