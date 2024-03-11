# Raman Spectra Analysis Tool

This project provides a program capable of fully analyzing Raman and IR spectra from different types of spectrometers. It features a graphical user interface (GUI) and supports various baseline removal and peak finding strategies.

## Features

- **Peak Fitting**: The tool allows for fitting of peaks with full Voigt profiles, pseudo Voigt profiles, Fano, Fano Voigt, Gaussian and Lorentzian models individually. It now incorporates asymmetric curve fitting, using Bimodal Gauss Gauss, Bimodal Lorentz-Lorentz, Mixed bimodal Lorentz Gauss, and more complex algorithms like Pearson IV and Sigmoidal Gauss Lorentz.
- **Post-Processing**: Provides post-processed data between all peaks in the spectrum.
- **Batch Processing**: Includes a batch processing routine for handling multiple datasets at once.
- **Spectrometer Support**: Currently supports Horiba and B&WTech spectrometers, also Brukker IR data can be proceed. Any data file in .txt or .dat format can be processed, but samples are needed to create the correct parsing. Please send a sample if your data is not coming from a supported spectrometer.

## Future Developments

The project is under active development. Future functionalities will include map processing and expansion to IR data, CL data, and others.

## Contributing

If you want support for a particular spectrometer model, please send a sample to the author so the parsing can be included in the software.

## License

This project is licensed under the terms of the GNU General Public License.

Enjoy using the Raman Spectra Analysis Tool!
