# RFI Masker
Tiny tool to apply a numpy boolean array to the spectral flags column of a measurement set to mask out known RFI sources like GSM and satelite in L-band. Typically flagging known RFI sources produces much better backgrounds and tools like the mighty AOFlagger RFI hammer needs fewer iterations to converge to good solutions.

# To develop
```
python setup.py develop --user
(assumes .local/bin is in your PATH and .local/lib/python2.7/site-packages is in your PYTHONPATH)
```
# Users
## Preferred run method
```
docker build -t rfi_masker (inside cloned folder)
docker run -it --rm rfi_masker (anywhere)
```
## Install (if you're yet to discover the pure awesomeness of Docker)
```
python setup.py install --user
(assumes .local/bin is in your PATH and .local/lib/python2.7/site-packages is in your PYTHONPATH)
```
