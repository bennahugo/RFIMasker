# RFI Masker
    Copyright (C) 2016  SKA South Africa

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

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
