from pylab import *
from astropy.cosmology import Planck15
import astropy.units as u
import astropy.constants as consts
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord

class HSC_toolkits:

    ####
    # Mac
    # filter_path='/Users/koki/Dropbox/Work/HSC-project/HSC-tomography/filters/'
    # fitsmask_path='/Users/koki/Dropbox/Work/HSC-project/data/CHORUS-PDR1/CHORUS-PDR1_ancillary/CHORUS-PDR1_fitsmask/'

    # ubuntu
    # fitsmask_path='/home/koki/HSC-tomography/data/CHORUS-PDR1/CHORUS-PDR1_ancillary/CHORUS-PDR1_fitsmask/'
    # filter_path='/home/koki/HSC-tomography/filters/'

    # MacBook Pro (DAWN)
    fitsmask_path='/Users/koki/Projects/HSC-tomography/data/CHORUS-PDR1/CHORUS-PDR1_ancillary/CHORUS-PDR1_fitsmask/'
    filter_path='/Users/koki/Projects/HSC-tomography/filters/'

    ####

    # load specz catalog
    def load_specz_catalog(self,filename='/home/koki/HSC-SSP/data/DR3_CATALOGS/pdr3_dud_specz+photometry_catalog.csv.gz'):
        # load spec-z catalogue associated with HSC-SSP DR3
        print('loading spec-z catalog:  ', filename)
        self.specz=pd.read_csv(filename)
        return self.specz
    # load sky objects
    def load_sky_catalog(self,filename='/home/koki/HSC-SSP/data/DR3_CATALOGS/pdr3_sky_objects_tract9813_masked.csv.gz'):
        # load sky object catalogue associated with HSC-SSP DR3
        print('loading sky object catalog:  ', filename)
        self.sky=pd.read_csv(filename)
        return self.sky
    # load image data
    def load_image(self,filter='NB0718'):
        if hasattr(self,'image') == False:
            self.image={}
        if filter=='HSC-G' or filter=='HSC-R' or filter=='HSC-I' or filter=='HSC-Z' or filter=='HSC-Y':
            path='/home/koki/HSC-SSP/data/DR3/'
            fname=filter+'-9813-mosaic_all.fits'
        if filter=='NB0816':
            path='/home/koki/HSC-SSP/data/DR3/'
            fname=filter+'-9813-mosaic_all.fits'
        if filter=='NB0921':
            path='/home/koki/HSC-SSP/data/DR3/'
            fname=filter+'-9813-mosaic_all.fits'
        if filter=='NB0718':
            path='/home/koki/HSC-SSP/data/CHORUS-PDR1/'
            fname=filter+'-9813-mosaic_all.fits'
        self.image[filter]=fits.open(path+fname)
        print('loaded fits image to self.image[\'filter\'] from:  ', path+fname)

    # load spec-z catalogue associated with HSC-SSP DR3
    def load_SILVERRUSH_S18A_CATALOG(self,path='/home/koki/HSC-SSP/data/catalog_SSP_S18A_LAE/',
                                     catalog_name='cat_sext_lae_NB0816_udcosmos_after_visual_check_v1_PatchOverlapRemoved.cat'):
        # define column names
        f=open(path+'column_sext.column','r')
        column_names=f.readlines()
        f.close()
        column_names=[name.replace('\n','') for name in column_names]
        # read catalog
        df=pd.read_csv(path+catalog_name, delim_whitespace=True, names=column_names)
        self.LAE_catalog=df
        return self.LAE_catalog

    # load GOLDRUSH catalog
    def load_GOLDRUSH_S18A_CATALOG(self,path='/home/koki/HSC-tomography/data/goldrush_ver20210826/',
                                   catalog_name='all_riz_UD_COSMOS.cat'):
        # define column names
        f=open(path+catalog_name,'r')
        column_names=[]
        for line in f:
            if line.startswith('#'):
                column_names.append(line)
        f.close()
        column_names=[name.replace('\n','') for name in column_names]
        df=pd.read_csv(path+catalog_name, comment='#',delim_whitespace=True, names=column_names)
        self.GOLDRUSH_catalog=df
        return self.GOLDRUSH_catalog

    def select_background_objects(self,zmin=5.00,zmax=5.85,filter='NB0718',
                                  magcut=25.0, magcut_flag=False, specz_flag=True):
        if hasattr(self,'specz') == False:
            self.load_specz_catalog()
        if hasattr(self,'image') == False:
            self.load_image(filter=filter)
        # choose redsift range
        in_redshift = ( (self.specz['specz_redshift']>=zmin)
                      & (self.specz['specz_redshift']<=zmax) )
        # choose region
        wcs=WCS(self.image[filter][0].header)
        NAXIS1=self.image[filter][0].header['NAXIS1']
        NAXIS2=self.image[filter][0].header['NAXIS2']
        fov_edge1=wcs.pixel_to_world(0,0)
        fov_edge2=wcs.pixel_to_world(NAXIS1,NAXIS2)
        in_region = ( (self.specz['specz_ra'] *u.deg >= fov_edge2.ra)  &
                      (self.specz['specz_dec']*u.deg >= fov_edge1.dec) &
                      (self.specz['specz_ra'] *u.deg <= fov_edge1.ra)  &
                      (self.specz['specz_dec']*u.deg <= fov_edge2.dec) )
        # select objects
        idx = in_redshift & in_region

        if specz_flag:
            is_secure = self.specz['specz_flag_homogeneous'] == True
            idx = idx & is_secure

        if magcut_flag:
            is_bright = self.specz['z_apertureflux_15_mag'] <= magcut
            idx = idx & is_bright

        return self.specz[idx]

    def load_deimos10k_catalog(self,filepath='/home/koki/HSC-SSP/data/DEIMOS10k/DEIMOS10k_full_catalog.csv'):
        deimos10k_catalog=pd.read_csv(filepath, delim_whitespace=True)
        self.deimos10k=deimos10k_catalog

    def get_deimos_spectrum(self,specz_name,spec_path='/home/koki/HSC-SSP/data/DEIMOS10k/'):
        if hasattr(self,'deimos10k') == False:
            self.load_deimos10k_catalog()
        deimos_id=specz_name.replace('DEIMOS_2018_','')
        idx = self.deimos10k['ID'] == deimos_id
        ascii1d=self.deimos10k['ascii1d'][idx]
        print(deimos_id,ascii1d)
        if self.deimos10k['ascii1d'][idx].notna().values[0]:
            ascii1d_filename=ascii1d.to_numpy()[0]
            filepath=ascii1d_filename.replace('/data/COSMOS/spectra/deimos/', spec_path)
            spec1d=pd.read_csv(filepath,skiprows=5,names=['wavelength','flux','ivar'],
                               delim_whitespace=True)
            return spec1d, ascii1d_filename.replace('/data/COSMOS/spectra/deimos/', '')
        else:
            return zeros(3), 'No Data'

    def get_vuds_spectrum(self,vuds_id,spec_path='/home/koki/HSC-SSP/data/VUDS/cesam_VUDS-DR1_spectra_2021-09-23T022609Z/VUDS-COSMOS-DR1/'):
        import glob
        filename =glob.glob(spec_path+'/spec1d/'+'*'+str(vuds_id)+'*.fits')
        if filename:
            print('Opening ... : ', filename[0])
            hdu=fits.open(filename[0])
            spec1d={}
            header=hdu[0].header
            flux=hdu[0].data
            wavelength=header['CRVAL1']+arange(header['NAXIS1'])*header['CDELT1']
            spec1d['flux']=flux
            spec1d['wavelength']=wavelength

        noisename=glob.glob(spec_path+'/spec1dnoise/'+'*'+str(vuds_id)+'*.fits')
        if noisename:
            hdu=fits.open(noisename[0])
            noise=hdu[0].data
            spec1d['noise']=noise[0]

        return spec1d, filename[0].replace(spec_path,'')

    def count2flux(self,count,zero_point=27.0):
        # HSC zero-point 27.0 mag/DN (ref: HSC-SSP data release FAQ)
        flux0=10**(-(zero_point+48.60)/2.5) * u.erg/u.s/u.cm**2/u.Hz
        flux=flux0*count
        return flux
    def flux2mag(self,flux):
        if flux.to('Jy').value > 0.0 :
            ABmag=-2.5*log10(flux.to('Jy').value/3631.)
        else:
            ABmag=np.nan
        return ABmag
    def mag2flux(self,mag):
        flux=10**(-(mag+48.60)/2.5) * u.erg/u.s/u.cm**2/u.Hz
        return flux

    def plot_image_cutout(self,image,sky_aperture,cutout=5.0*u.arcsec,
                          plot_bkg_aperture=False,bkg_annulus=None):
        from astropy.nddata import Cutout2D
        # aperture position
        wcs=WCS(image.header)
        aperture=sky_aperture.to_pixel(wcs)
        x,y=aperture.positions
        # cutout size
        pixel_scale = image.header['CD2_2']*u.deg
        cutsize=(cutout/pixel_scale).to('').value
        cutout=Cutout2D(image.data,(x,y),size=(cutsize,cutsize))
        # plot
        min_count=-0.085
        max_count=+0.155
        fig,ax=plt.subplots(figsize=(4,3))
        im=ax.imshow(cutout.data,origin='lower')
        fig.colorbar(im)
        circ=plt.Circle((cutsize/2,cutsize/2),aperture.r,edgecolor='white',facecolor='none')
        ax.add_patch(circ)
        if plot_bkg_aperture:
            bkg_aperture=bkg_annulus.to_pixel(wcs)
            circ1=plt.Circle((cutsize/2,cutsize/2),bkg_aperture.r_in, edgecolor='red',facecolor='none')
            circ2=plt.Circle((cutsize/2,cutsize/2),bkg_aperture.r_out,edgecolor='red',facecolor='none')
            ax.add_patch(circ1)
            ax.add_patch(circ2)
        plt.show()

    def measure_photometry(self,ra,dec,tract='9813',patch='0,0',filter='HSC-Z',
                           diameter=2.0*u.arcsec,path='/media/koki/MAC_DATA/HSC-SSP/',
                           use_mask=True,plot_cutout=False,
                           local_bkg_subt=False,r_in=2.0*u.arcsec,r_out=4.0*u.arcsec):
        from astropy.coordinates import SkyCoord
        from photutils.aperture import SkyCircularAperture, SkyCircularAnnulus
        from photutils.aperture import CircularAperture, CircularAnnulus,aperture_photometry, ApertureStats
        from astropy.stats import SigmaClip
        # check if data exists
        if tract=='None' or patch=='None':
            print('No tract or patch exists. Skipping photometry (=np.nan) ')
            flux=np.nan * u.erg/u.s/u.cm**2/u.Hz
            flux_err=np.nan * u.erg/u.s/u.cm**2/u.Hz
            ABmag=nan
            return (flux, flux_err, ABmag)
        # open fits file
        image_file='calexp-'+filter+'-'+tract+'-'+patch+'.fits'
        if (filter=='HSC-G' or filter=='HSC-R' or filter=='HSC-I' or
            filter=='HSC-Z' or filter=='HSC-Y' or filter=='NB0387' or
            filter=='NB0816' or filter=='NB0921' or filter=='NB1010'):
            hdul=fits.open(path+'/DR3/'+image_file)
        elif (filter=='NB0718'):
            hdul=fits.open(path+'/CHORUS-PDR1/'+image_file)
            if use_mask:
                mask_path='/home/koki/HSC-SSP/data/CHORUS-PDR1/CHORUS-PDR1_ancillary/CHORUS-PDR1_fitsmask/CombMask_NB0718_GAIA_BRIGHT_STARS/'
                mask_file='CombMask_NB0718_GAIA_BRIGHT_STARS_t'+tract+'p'+patch+'.fits'
                fitsmask=fits.open(mask_path+mask_file)
                mask=(fitsmask[1].data == 1) # 1 (=True) for masked region, 0 (=False) for unmasked region
        else:
            print('HSC data not found! You need to download the data!')
        # obtain IMAGE, MASK, VARINCE from HSC Coadds
        image=hdul[1] # image
        #mask=hdul[2]  # mask
        var=hdul[3]   # variance
        # define sky circular aperture and do photometry
        sky_position=SkyCoord(ra,dec,unit='deg',frame='icrs')
        sky_aperture=SkyCircularAperture(sky_position, r=diameter/2)
        # define sky circular annulus to estimate local background
        if use_mask & (filter=='NB0718'):
            print('photometry with sky mask regions for NB0718...')
            if local_bkg_subt: # photometry with local background subtraction
                bkg_annulus=SkyCircularAnnulus(sky_position, r_in=r_in, r_out=r_out)
                sigclip=SigmaClip(sigma=3.0,maxiters=100)
                local_bkg=ApertureStats(image.data, bkg_annulus, mask=mask, wcs=WCS(image.header), sigma_clip=sigclip)
                # photometry of an object - local background level
                phot_table=aperture_photometry(image.data-local_bkg.median, sky_aperture, error=sqrt(var.data), mask=mask, wcs=WCS(image.header))
                phot_flag =aperture_photometry(fitsmask[1].data, sky_aperture, wcs=WCS(image.header))
            else: # photometry without local background subtraction
                # photometry of an object
                phot_table=aperture_photometry(image.data, sky_aperture, error=sqrt(var.data), mask=mask, wcs=WCS(image.header))
                phot_flag =aperture_photometry(fitsmask[1].data, sky_aperture, wcs=WCS(image.header))
            if phot_flag['aperture_sum'][0]>=1:
                self.nb718_mask_flag=True
            else:
                self.nb718_mask_flag=False
        else:
            if local_bkg_subt: # photometry with local background subtraction
                bkg_annulus=SkyCircularAnnulus(sky_position, r_in=r_in, r_out=r_out)
                sigclip=SigmaClip(sigma=3.0,maxiters=100)
                local_bkg=ApertureStats(image.data, bkg_annulus, wcs=WCS(image.header), sigma_clip=sigclip)
                # photometry of an object - local background level
                phot_table=aperture_photometry(image.data-local_bkg.median, sky_aperture, error=sqrt(var.data), wcs=WCS(image.header))
            else:
                # photometry of an object
                phot_table=aperture_photometry(image.data, sky_aperture, error=sqrt(var.data), wcs=WCS(image.header))
        # save results
        count     = phot_table['aperture_sum'][0]
        count_err = phot_table['aperture_sum_err'][0]
        flux      = self.count2flux(count)
        flux_err  = self.count2flux(count_err)
        ABmag     = self.flux2mag(flux)
        if plot_cutout:
            if local_bkg_subt:
                self.plot_image_cutout(image,sky_aperture,cutout=30.*u.arcsec,
                                plot_bkg_aperture=True,bkg_annulus=bkg_annulus)
            else:
                self.plot_image_cutout(image,sky_aperture,cutout=30.*u.arcsec)
        return (flux, flux_err, ABmag)

    def find_tract_patch(self,ra,dec,unit='deg',frame='icrs',verbose=True,
                         path='/media/koki/MAC_DATA/HSC-SSP/DR3/'):
        from glob import glob
        object=SkyCoord(ra,dec,unit=unit,frame=frame)
        # HSC-Z: reference filter image to get all tract and patch number
        all_files=glob(path+'calexp-HSC-Z-*.fits')
        notract=True # flag to check if ra,dec are found in downloaded images
        for file in all_files:
            header=fits.getheader(file,1) # get header of image HDU
            w=WCS(header)
            if w.footprint_contains(object): # check if the object is inside the footprint of the image
                tract_patch_xy=file.replace(path+'calexp-HSC-Z-','')  # get tract and patch number
                tract_patch_xy=tract_patch_xy.replace('.fits','')
                tract, patch_s = tract_patch_xy.split(sep='-')
                notract = False
                break
        if notract:
            print('No tract found. Need to download data if needed')
            tract='None'
            patch_s='None'
        if verbose:
            print('replacing to (tract, patch_s) = ', '(', tract,',',patch_s,')')
        return (tract,patch_s)

    # Subaru/HSC property, filter FWHM etc
    def filter(self,filter='NB718',filter_path=filter_path,#'/Users/koki/Dropbox/Work/HSC-project/HSC-tomography/filters/',
                    include_throughput=True,verbose=True):
        from scipy.signal import find_peaks, peak_widths
        from scipy.interpolate import interp1d
        filter_name='wHSC-'+filter+'.txt'
        if verbose: print('Getting the filter properties from ... ', filter_path+filter_name)
        filter={}
        wavelength,filter_transmission=genfromtxt(filter_path+filter_name,unpack=True)
        filter['wavelength']=wavelength
        filter['filter_transmission']=filter_transmission
        if include_throughput:
            wl,ccd_efficiency=genfromtxt(filter_path+'qe_ccd_HSC.txt',unpack=True)
            ccd_efficiency=interp1d(wl,ccd_efficiency,bounds_error=False,fill_value=0.0)
            wl,dewar_window_transmission=genfromtxt(filter_path+'throughput_win.txt',unpack=True)
            dewar_window_transmission=interp1d(wl,dewar_window_transmission,bounds_error=False,fill_value=0.0)
            wl,primary_focus_transmission=genfromtxt(filter_path+'throughput_popt2.txt',unpack=True)
            primary_focus_transmission=interp1d(wl,primary_focus_transmission,bounds_error=False,fill_value=0.0)

            filter_transmission=filter_transmission*ccd_efficiency(wavelength)*dewar_window_transmission(wavelength)*primary_focus_transmission(wavelength)
            filter['filter_transmission']=filter_transmission
        z_lya=wavelength/1215.67-1
        z_lya=interp1d(arange(z_lya.size),z_lya)
        wavelength=interp1d(arange(wavelength.size),wavelength)
        peak=where(filter_transmission==filter_transmission.max())[0]
        widths,width_heights,left_ips,right_ips=peak_widths(filter_transmission, peak, rel_height=0.5)
        z_min=z_lya(left_ips)[0]
        z_max=z_lya(right_ips)[0]
        wavelength_peak=wavelength(peak)[0]
        wavelength_min =wavelength(left_ips)[0]
        wavelength_max =wavelength(right_ips)[0]
        wavelength_cen =(wavelength_max+wavelength_min)/2
        filter['z_peak']=wavelength_peak/1215.67-1
        filter['z_cen'] =wavelength_cen/1215.67-1
        filter['z_min']=z_min
        filter['z_max']=z_max
        filter['dz']=z_max-z_min
        filter['wavelength_min']=wavelength_min
        filter['wavelength_max']=wavelength_max
        filter['wavelength_peak']=wavelength_peak
        filter['wavelength_cen']=wavelength_cen
        filter['wavelength_FWHM']=wavelength_max-wavelength_min
        return filter
    
    # JWST NIRCam filter transmission
    def NIRCam_filter(self,filter='F115W',filter_path='/Users/koki/Projects/HSC-tomography/filters/nircam_throughput/',
                           output='mean_throughputs',verbose=True):
        from scipy.interpolate import interp1d
        from astropy.table import Table
        filter_name=filter+'_mean_system_throughput.txt'
        filename=filter_path+'/'+output+'/'+filter_name
        if verbose: print('Getting the filter properties from ... ',filename)
        filter={}
        dat=Table.read(filename,format='ascii')
        filter['wavelength']=dat['Microns']*u.um
        filter['filter_transmission']=dat['Throughput']
        return filter

    # UltraVISTA YJHKs filter transmission
    def VISTA_filter(self,filter='Ks',filter_path='/Users/koki/Projects/HSC-tomography/filters/VISTA_throughput/',
                           verbose=True):
        from scipy.interpolate import interp1d
        filter_name='VISTA_Filters_at80K_forETC_'+filter+'.dat'
        filename=filter_path+'/'+filter_name
        if verbose: print('Getting the filter properties from ... ',filename)
        filter={}
        wavelength,filter_transmission=genfromtxt(filename,unpack=True)
        filter['wavelength']=wavelength*u.nm
        filter['filter_transmission']=filter_transmission/100
        return filter

    # Euclid YJH filter transmission
    def Euclid_filter(self,filter='Y',filter_path='/Users/koki/Projects/HSC-tomography/filters/Euclid_throughput/',
                           verbose=True):
        from scipy.interpolate import interp1d
        filter_name='Euclid-'+filter+'.dat'
        filename=filter_path+'/'+filter_name
        if verbose: print('Getting the filter properties from ... ',filename)
        filter={}
        dat=Table.read(filename,format='ascii')
        filter['wavelength']=dat['WAVE']*u.nm
        filter['filter_transmission']=dat['T_TOTAL']
        return filter

    def mag2UVmag(self,BBmag,redshift,flag99=True):
        from astropy.cosmology import Planck15
        Kcorrection=2.5*np.log10(1+redshift)
        Muv=BBmag-Planck15.distmod(redshift).to('mag').value+Kcorrection
        if flag99: Muv[(BBmag==99.)]=np.nan
        return Muv
    def mag2lyaflux(self,NBmag, BBmag, filtername='NB718',include_throughput=True):
        from scipy.interpolate import interp1d
        filter=self.filter(filter=filtername,include_throughput=include_throughput,verbose=False)
        NBflux=self.mag2flux(NBmag)
        BBflux=self.mag2flux(BBmag)
        nu=consts.c/(filter['wavelength']*u.angstrom)
        factor=filter['filter_transmission'].max()/abs(trapz(filter['filter_transmission'],x=nu))
        lyaflux=(NBflux-BBflux)/factor
        return lyaflux.to('erg/(s cm2)')
    def lyaflux2luminosity(self,lyaflux,filtername='NB718',include_throughput=True):
        from astropy.cosmology import Planck15
        filter=self.filter(filter=filtername,include_throughput=include_throughput,verbose=False)
        z=filter['z_peak']
        DL=Planck15.luminosity_distance(z)
        La=lyaflux * (4*pi*DL**2)
        return (La.to('erg/s'), z)

    def get_mask(self,xg,yg,fname='CombMask_NB0718-mosaic_all_GAIA_BRIGHT_STARS.fits',
                 path=fitsmask_path):
                 #'/home/koki/HSC-SSP/data/CHORUS-PDR1/CHORUS-PDR1_ancillary/CHORUS-PDR1_fitsmask/'):
        " Get resampled mask region at (RA,DEC)=xg,yg "
        from scipy.interpolate import interpolate
        if hasattr(self,'mask') == False:
            hdu=fits.open(path+fname)
            header=hdu[0].header
            dat=hdu[0].data    # masked region data: =1 if masked, =0 if not masked
            dat=1-dat          # swap 0 -> 1 and 0 -> 1
            wcs=WCS(header)
            pix1=np.arange(header['NAXIS1'])
            pix2=np.arange(header['NAXIS2'])
            coord=wcs.pixel_to_world(pix1,pix2)
            x=coord.ra.value
            y=coord.dec.value
            pix_area=np.abs(header['CDELT1']*u.Unit(header['CUNIT1'])*header['CDELT2']*u.Unit(header['CUNIT2']))
            self.survey_area=np.sum(dat)*pix_area
            self.mask=interpolate.interp2d(x,y,dat,kind='linear',bounds_error=False,fill_value=0)
        return self.mask(xg,yg)

    def make_random_catalog(self,size=10000,fname='CombMask_NB0718-mosaic_all_GAIA_BRIGHT_STARS.fits',
                            path=fitsmask_path):
#                            path='/home/koki/HSC-SSP/data/CHORUS-PDR1/CHORUS-PDR1_ancillary/CHORUS-PDR1_fitsmask/'):
        print('making random objects in the footprint of ...', fname)
        hdu=fits.open(path+fname)
        header=hdu[0].header
        wcs=WCS(header)
        dat=hdu[0].data    # masked region data: =1 if masked, =0 if not masked
        idx_x,idx_y=np.where(dat==0)
        idx_length=np.arange(len(idx_x))
        samples=np.random.choice(idx_length,size=size,replace=True)
        random_objects=wcs.array_index_to_world(idx_x[samples],idx_y[samples])
        return random_objects

    def load_MCMC_result(self,
            source='LAE', # or DEIMOS10k
            catalog='NB0718_bkg_LAE_catalog_SSP_S18A_ALL_HSC-PHOTOMETRY-LIMMAG_WITH_MASK_Z5SIGMA.csv',
            MCMC_directory='MCMC_results_SSP-S18A/',
            apply_mask_flag=True,
            apply_g_cut=False, g_sigma=3.0,
            apply_nb718_z_color_cut=False,
            apply_upper_outlier_cut=False,Pr_threshold=0.68,
            limit_redshift=False, z_min=4.98, z_max=5.89):
        import glob
        # get catalog
        objects=pd.read_csv('catalogs/'+catalog)
        if source=='LAE':
            objects=objects.rename(columns={'#2 [NB]  NUMBER Running object number': 'name'})
        if source=='DEIMOS10k':
            objects=objects.rename(columns={'specz_name': 'name'})
        if apply_mask_flag:
            idx=(objects['nb718_mask_flag15']==False)
            objects=objects[idx]
        if apply_g_cut: # require non-detection in g
            g_cut=(objects['g_flux15']<g_sigma*objects['g_limflux15'])
            objects=objects[g_cut]
        if apply_nb718_z_color_cut: # require NB718-z color must be red
            nb718_cut=objects['nb718_flux15']-objects['z_flux15']<0.0
            objects=objects[nb718_cut]
        if limit_redshift:
            redshift_range=(objects['specz_redshift']>z_min)&(objects['specz_redshift']<z_max)
            objects=objects[redshift_range]
        # get all computed MCMC chains
        files=glob.glob(MCMC_directory+"/*.csv")
        All_MCMC_Chains={}
        for file in files:
            name=file.replace(MCMC_directory,'').replace('_MCMC_result.csv','')
            All_MCMC_Chains[name]=pd.read_csv(file)
        # store reduced results ONLY for the selected objects in the catalog
        Result={'name':[], 'RA':[], 'DEC':[],
                'TIGM_mean': [], 'TIGM_var': [],'TIGM_median': [],
                'TIGM_05tile': [], 'TIGM_16tile': [], 'TIGM_84tile': [], 'TIGM_95tile': [],
                'beta_mean': [], 'beta_median': [], 'beta_var': [],
                'beta_05tile': [], 'beta_16tile': [], 'beta_84tile': [],'beta_95tile': [],
                'Muv_mean': [], 'Muv_median': [], 'Muv_var': [],
                'Muv_05tile': [], 'Muv_16tile': [], 'Muv_84tile': [],'Muv_95tile': []
                }
        for index,obj in objects.iterrows():
            name,ra,dec=obj[['name','ra','dec']]
            name=str(name)
            if name in All_MCMC_Chains.keys():
                # TIGM
                TIGM_mean=All_MCMC_Chains[name]['TIGM'].mean()
                TIGM_var =All_MCMC_Chains[name]['TIGM'].var()
                TIGM_median=np.median(All_MCMC_Chains[name]['TIGM'])
                TIGM_05tile=np.quantile(All_MCMC_Chains[name]['TIGM'],0.05)
                TIGM_16tile=np.quantile(All_MCMC_Chains[name]['TIGM'],0.16)
                TIGM_84tile=np.quantile(All_MCMC_Chains[name]['TIGM'],0.84)
                TIGM_95tile=np.quantile(All_MCMC_Chains[name]['TIGM'],0.95)
                # beta
                beta_mean=All_MCMC_Chains[name]['beta'].mean()
                beta_median=np.median(All_MCMC_Chains[name]['beta'])
                beta_var =All_MCMC_Chains[name]['beta'].var()
                beta_05tile=np.quantile(All_MCMC_Chains[name]['beta'],0.05)
                beta_16tile=np.quantile(All_MCMC_Chains[name]['beta'],0.16)
                beta_84tile=np.quantile(All_MCMC_Chains[name]['beta'],0.84)
                beta_95tile=np.quantile(All_MCMC_Chains[name]['beta'],0.95)
                # Muv
                Muv_mean=All_MCMC_Chains[name]['Muv'].mean()
                Muv_median=np.median(All_MCMC_Chains[name]['Muv'])
                Muv_var =All_MCMC_Chains[name]['Muv'].var()
                Muv_05tile=np.quantile(All_MCMC_Chains[name]['Muv'],0.05)
                Muv_16tile=np.quantile(All_MCMC_Chains[name]['Muv'],0.16)
                Muv_84tile=np.quantile(All_MCMC_Chains[name]['Muv'],0.84)
                Muv_95tile=np.quantile(All_MCMC_Chains[name]['Muv'],0.95)
                # store the result
                if apply_upper_outlier_cut:
                    upper_outliers=All_MCMC_Chains[name]['TIGM']>1.0
                    #lower_outliers=All_MCMC_Chains[name]['TIGM']<0.0
                    Pr_upper=len(All_MCMC_Chains[name]['TIGM'][upper_outliers])/len(All_MCMC_Chains[name]['TIGM'])
                    #Pr_lower=len(All_MCMC_Chains[name]['TIGM'][lower_outliers])/len(All_MCMC_Chains[name]['TIGM'])
                    if (Pr_upper>Pr_threshold):# or (Pr_lower>Pr_threshold):
                        print('reject outlier :  '+name+'  Pr(TIGM>1)=%.3f' % Pr_upper)
                    else:
                        Result['name'].append(name)
                        Result['RA'].append(ra)
                        Result['DEC'].append(dec)
                        Result['TIGM_mean'].append(TIGM_mean)
                        Result['TIGM_median'].append(TIGM_median)
                        Result['TIGM_05tile'].append(TIGM_05tile)
                        Result['TIGM_16tile'].append(TIGM_16tile)
                        Result['TIGM_84tile'].append(TIGM_84tile)
                        Result['TIGM_95tile'].append(TIGM_95tile)
                        Result['TIGM_var'].append(TIGM_var)
                        Result['beta_mean'].append(beta_mean)
                        Result['beta_median'].append(beta_median)
                        Result['beta_var'].append(beta_var)
                        Result['beta_05tile'].append(beta_05tile)
                        Result['beta_16tile'].append(beta_16tile)
                        Result['beta_84tile'].append(beta_84tile)
                        Result['beta_95tile'].append(beta_95tile)
                        Result['Muv_mean'].append(Muv_mean)
                        Result['Muv_median'].append(Muv_median)
                        Result['Muv_var'].append(Muv_var)
                        Result['Muv_05tile'].append(Muv_05tile)
                        Result['Muv_16tile'].append(Muv_16tile)
                        Result['Muv_84tile'].append(Muv_84tile)
                        Result['Muv_95tile'].append(Muv_95tile)
                else:
                    Result['name'].append(name)
                    Result['RA'].append(ra)
                    Result['DEC'].append(dec)
                    Result['TIGM_mean'].append(TIGM_mean)
                    Result['TIGM_median'].append(TIGM_median)
                    Result['TIGM_05tile'].append(TIGM_05tile)
                    Result['TIGM_16tile'].append(TIGM_16tile)
                    Result['TIGM_84tile'].append(TIGM_84tile)
                    Result['TIGM_95tile'].append(TIGM_95tile)
                    Result['TIGM_var'].append(TIGM_var)
                    Result['beta_mean'].append(beta_mean)
                    Result['beta_median'].append(beta_median)
                    Result['beta_var'].append(beta_var)
                    Result['beta_05tile'].append(beta_05tile)
                    Result['beta_16tile'].append(beta_16tile)
                    Result['beta_84tile'].append(beta_84tile)
                    Result['beta_95tile'].append(beta_95tile)
                    Result['Muv_mean'].append(Muv_mean)
                    Result['Muv_median'].append(Muv_median)
                    Result['Muv_var'].append(Muv_var)
                    Result['Muv_05tile'].append(Muv_05tile)
                    Result['Muv_16tile'].append(Muv_16tile)
                    Result['Muv_84tile'].append(Muv_84tile)
                    Result['Muv_95tile'].append(Muv_95tile)
        # convert all lists to arrays
        for key in Result:
            Result[key]=np.array(Result[key])
        return (Result, All_MCMC_Chains)


    def tau_eff(self,z,type='Schmidt19'):
        if type=='Schmidt19':
            tau_eff=0.00126*exp(3.294*sqrt(z))
        if type=='Becker13':
            tau0=0.751
            beta=2.90
            C=-0.132
            z0=3.5
            tau_eff=tau0*((1+z)/(1+z0))**beta+C
        if type=='Bosman21':
            tau0=0.30
            beta=13.7
            C=1.35
            z0=4.8
            tau_eff=tau0*((1+z)/(1+z0))**beta+C
        return tau_eff
    # set angular bins
    def angular_bins(self,min_bin,max_bin,bins,type='linear-bin'):
        #theta_bins,dtheta_bin=np.linspace(min_bin,max_bin,bins,retstep=True)
        #return theta_bins, dtheta_bin
        if type=='linear-bin':
            theta_edges=np.linspace(min_bin,max_bin,bins+1)
        if type=='log-bin':
            theta_edges=np.logspace(np.log10(min_bin.to('arcmin').value),
                                    np.log10(max_bin.to('arcmin').value),
                                    bins+1)
            theta_edges=theta_edges*u.arcmin
        theta_bins=np.zeros(bins)*u.arcmin
        for i in range(bins):
            theta_bins[i]=(theta_edges[i]+theta_edges[i+1])/2
        return theta_bins, theta_edges

    # estimator: mean Lya forest transmission
    def mean_transmission(self, quantity=None):
        return quantity.mean()

    # estimator: mean Lya transmitted flux profile around galaxies
    def mean_profile_around_objects(self,
                            foreground_objects, background_objects, weights=None,
                            min_bin=1*u.arcmin,max_bin=30*u.arcmin,bins=10,
                            type='linear-bin'):
        """
        Input: foreground objects : coord of fg. objects (astropy.SkyCoord)
               background_objects : coord of bg. objects (astropy.SkyCoord)
               weights            : quantity to be averaged (numpy.array)
        """
        # angular bins
        theta_bins, theta_edges=self.angular_bins(min_bin,max_bin,bins,type=type)
        # define arrays for mean flux profile and number of pairs
        mean_flux=np.zeros(bins)
        npairs=np.zeros(bins)
        # compute mean TIGM around foreground_objects
        print('computing an angular mean profile ...')
        for foreground_object in foreground_objects:
            theta=foreground_object.separation(background_objects)
            for n in range(bins):
                in_bin=( (theta>=theta_edges[n]) & (theta< theta_edges[n+1]) )
                mean_flux[n]+=np.sum(weights[in_bin])
                npairs[n]+=weights[in_bin].size
        mean_flux=mean_flux/npairs
        return(mean_flux, npairs, theta_bins)
    # max posterior estimator: mean Lya transmitted flux profile around galaxies
    def maximum_posterior_estimation_of_mean_profile_around_objects(self,
                    foreground_objects, background_objects,
                    object_names=None, MCMC_chains=None, random_size=100,
                    min_bin=1*u.arcmin,max_bin=30*u.arcmin,bins=10,
                    type='linear-bin'):
        """
        Input: foreground objects : coord of fg. objects (astropy.SkyCoord)
               background_objects : coord of bg. objects (astropy.SkyCoord)
               object_names:      : names of background bojects (numpy.str)
               MCMC_chains        : MCMC chains ('name','TIGM','Muv','beta')
                                    results from the IGM toolkits (numpy.dict)
        """
        print('Maximum posterior estimation of mean profile around objects ...')
        # angular bins
        theta_bins, theta_edges=self.angular_bins(min_bin,max_bin,bins,type=type)
        # define arrays for mean flux profile and number of pairs for each realization
        mean_flux=np.zeros([bins,random_size])
        npairs=np.zeros([bins,random_size])
        # compute mean TIGM around foreground_objects
        for foreground_object in foreground_objects:
            # generate weights by random sampling from individual posteriors
            weights=np.zeros([len(background_objects),random_size])
            for n, name in zip( range(len(background_objects)), object_names ) :
                weights[n,:]=np.random.choice( MCMC_chains[name]['TIGM'], size=random_size )
            # compute the profile by binning weights into angular bins
            theta=foreground_object.separation(background_objects)
            for n in range(bins):
                in_bin=( (theta>=theta_edges[n]) & (theta< theta_edges[n+1]) )
                for i in range(random_size):
                    mean_flux[n,i]+=np.sum(weights[in_bin,i])
                    npairs[n,i]+=weights[in_bin,i].size
        mean_flux=mean_flux/npairs
        return (mean_flux, npairs, theta_bins)

    # variance of the estimator: mean Lya transmitted flux profile around galaxies
    def variance_profile_around_objects(self,
                    foreground_objects, background_objects, weights=None,
                    min_bin=1*u.arcmin,max_bin=30*u.arcmin,bins=10,
                    type='linear-bin'):
        """
        Input: foreground objects : coord of fg. objects (astropy.SkyCoord)
               background_objects : coord of bg. objects (astropy.SkyCoord)
               weights            : quantity to be averaged (numpy.array)
        """
        # angular bins
        theta_bins, theta_edges=self.angular_bins(min_bin,max_bin,bins,type=type)
        variance=np.zeros(bins)
        npairs=np.zeros(bins)
        # compute mean TIGM around foreground_objects
        print('computing an angular variance profile ...')
        for foreground_object in foreground_objects:
            theta=foreground_object.separation(background_objects)
            for n in range(bins):
                in_bin=( (theta>=theta_edges[n]) & (theta< theta_edges[n+1]) )
                variance[n]+=np.sum(weights[in_bin])
                npairs[n]+=weights[in_bin].size
        variance=variance/(npairs**2)
        return(variance, npairs, theta_bins)
    # covariance of the estimator (analytic): mean Lya transmitted flux profile around galaxies
    def covariance_profile_around_objects(self,
                    foreground_objects, background_objects, weights=None,
                    min_bin=1*u.arcmin,max_bin=30*u.arcmin,bins=10,
                    type='linear-bin'):
        # angular bins
        theta1_bins, theta1_edges=self.angular_bins(min_bin,max_bin,bins,type=type)
        theta2_bins, theta2_edges=self.angular_bins(min_bin,max_bin,bins,type=type)
        covariance=np.zeros([bins,bins])
        dummy=np.zeros([bins,bins])
        npairs=np.zeros(bins)
        N1N2=np.zeros([bins,bins])
        for i in range(background_objects.size):
            theta=background_objects[i].separation(foreground_objects)
            for n in range(bins):
                in_bin1=( (theta>=theta1_edges[n]) & (theta< theta1_edges[n+1]) )
                npairs[n]+=foreground_objects[in_bin1].size
        # KK: I'm still not convinced by this analytic formula/calculation!! DEBUG/CHECK!!
        print('computing an angular covariance profile ...')
        for i in range(background_objects.size):
            theta=background_objects[i].separation(foreground_objects)
            for n in range(bins):
                for m in range(bins):
                    in_bin1=( (theta>=theta1_edges[n]) & (theta< theta1_edges[n+1]) )
                    in_bin2=( (theta>=theta2_edges[m]) & (theta< theta2_edges[m+1]) )
                    N1=foreground_objects[in_bin1].size
                    N2=foreground_objects[in_bin2].size
                    if n==m:
                        N1N2[n,m]=N1
                    else:
                        N1N2[n,m]=N1*N2
                    covariance[n,m]=covariance[n,m]+N1N2[n,m]*weights[i]
        for n in range(bins):
            for m in range(bins):
                covariance[n,m]=covariance[n,m]/(npairs[n]*npairs[m])

        return (covariance, npairs, npairs, theta1_bins, theta2_bins)

    # Jackknife error estimator: mean Lya transmitted flux profile around galaxies
    def jackknife_error_estimate(self,
                    foreground_objects, background_objects, weights=None,
                    random_objects=None, internal_random=True,
                    min_bin=1*u.arcmin,max_bin=30*u.arcmin,bins=10,
                    jk_region_size=20,plot_check=True,type='linear-bin',
                    return_covariance_matrix=False):
        print('Jackknife error estimation (cross-correlation function)...')
        # initialize random seed
        np.random.seed(seed=0)
        # define jackknife regions from k-mean clustering of random objects
        #  https://github.com/esheldon/kmeans_radec
        #  Note: python setup.py install is deprecated. 
        #   Use cd kmeans_radec
        #       pip install . 
        from kmeans_radec import KMeans, kmeans_sample 

        if internal_random: random_objects=self.make_random_catalog(size=50000)
        R=np.vstack( [random_objects.ra.value, random_objects.dec.value] ).T
        km=kmeans_sample(R, jk_region_size, maxiter=100, tol=1.0e-5)
        if plot_check:
            fig,ax=plt.subplots(figsize=(6,5))
            labels=km.find_nearest(R)
            ax.scatter(R[:,0],R[:,1],c=labels,s=1,marker='.',cmap='tab20')
            ax.set_aspect('equal')
            ax.set_xlabel('RA [deg]')
            ax.set_ylabel('DEC [deg]')
            ax.invert_xaxis()
            plt.tight_layout()
            plt.pause(0.05) #plt.show()
        # get the Jackknife labels for all objects
        jk_regions={}
        jk_regions['foreground_objects']=km.find_nearest(
                                            np.vstack( [ foreground_objects.ra.value,
                                                         foreground_objects.dec.value ] ).T
                                                        )
        jk_regions['background_objects']=km.find_nearest(
                                            np.vstack( [ background_objects.ra.value,
                                                         background_objects.dec.value ] ).T
                                                        )
        # angular bins
        theta_bins, theta_edges=self.angular_bins(min_bin,max_bin,bins,type=type)
        # define mean transmitted flux profiles for jackknifed sample
        jk_mean_flux=[]
        for k in range(jk_region_size):
            print('At ',k,'th Jackknife sample ...')
            # jackknife the foreground objects
            jackknife_sample_fg = (jk_regions['foreground_objects'] != k)
            jk_foreground_objects=foreground_objects[jackknife_sample_fg]
            # jackknife the background objects
            jackknife_sample_bg = (jk_regions['background_objects'] != k)
            jk_background_objects=background_objects[jackknife_sample_bg]
            jk_weights=weights[jackknife_sample_bg]
            # compute mean transmitted profile for the jackknifed sample
            mean_flux_k, _ , _ = self.mean_profile_around_objects(
                    jk_foreground_objects, jk_background_objects, weights=jk_weights,
                    min_bin=min_bin,max_bin=max_bin,bins=bins,type=type
                    )
            jk_mean_flux.append(mean_flux_k)
        # calculate jackknife error
        # np.var(x) is the average of the squared deviations from the mean,
        #   i.e., var = mean(x), where x = abs(a - a.mean())**2.
        jk_mean_flux=np.array(jk_mean_flux) # convert to array for easy computation
        var=np.var(jk_mean_flux,axis=0)  # = (1/N_JK) * Sum( [f(x)^k - <f(x)^k>]^2 )
        jackknife_variance=(jk_region_size-1)*var # = ((N_JK-1)/N_JK) * Sum( [f(x)^k - <f(x)^k>]^2 )

        if return_covariance_matrix:
            covariance=self.covariance_matrix(jk_mean_flux,type='Jackknife')
            return (jackknife_variance, theta_bins, covariance)
        else:
            return (jackknife_variance, theta_bins)

    # Bootstrapping error estimator: mean Lya transmitted flux profile around galaxies
    def bootstrap_error_estimate(self,
                    foreground_objects, background_objects, weights=None,
                    random_objects=None, internal_random=True, return_covariance_matrix=False,
                    min_bin=1*u.arcmin,max_bin=30*u.arcmin,bins=10,
                    subregion_size=20,resampling_size=20,plot_check=True,type='linear-bin'):
        print('Bootstrapping error estimation (cross-correlation function)...')
        # initialize random seed
        np.random.seed(seed=0)
        # define subvolume regions from k-mean clustering of random objects
        from kmeans_radec import KMeans, kmeans_sample
        if internal_random: random_objects=self.make_random_catalog(size=50000)
        R=np.vstack( [random_objects.ra.value, random_objects.dec.value] ).T
        km=kmeans_sample(R, subregion_size, maxiter=100, tol=1.0e-5, verbose=0)
        if plot_check:
            fig,ax=plt.subplots(figsize=(6,5))
            labels=km.find_nearest(R)
            ax.scatter(R[:,0],R[:,1],c=labels,s=1,marker='.')
            ax.set_aspect('equal')
            ax.set_xlabel('RA [deg]')
            ax.set_ylabel('DEC [deg]')
            plt.tight_layout()
            plt.pause(0.05) #plt.show()
        # get the subregion labels for all objects
        subregion_labels={}
        subregion_labels['foreground_objects']=km.find_nearest(
                                                np.vstack( [ foreground_objects.ra.value,
                                                             foreground_objects.dec.value ] ).T
                                                             )
        subregion_labels['background_objects']=km.find_nearest(
                                                np.vstack( [ background_objects.ra.value,
                                                             background_objects.dec.value ] ).T
                                                             )
        # bootstrap resampling the subregions
        resampled_regions=resampling_size*[None]
        for k in range(resampling_size):
            resampled_regions[k]=np.random.choice(
                                            np.unique(km.labels),
                                            size=subregion_size,
                                            replace=True
                                            )
        # angular bins
        theta_bins, theta_edges=self.angular_bins(min_bin,max_bin,bins,type=type)
        # define mean transmitted flux profiles for resampled regions
        resampled_mean_flux=[]
        for k in range(resampling_size):
            print('At ',k,'th Bootstrapping sample ...')
            # bootstrap resampling of the foreground objects
            bootstrap_sample=np.array([
                            (subregion_label in resampled_regions[k])
                            for subregion_label in subregion_labels['foreground_objects']
                            ])
            bootstrap_foreground_objects=foreground_objects[bootstrap_sample]
            # booststrap resampling of the background objects
            bootstrap_sample=np.array([
                            (subregion_label in resampled_regions[k])
                            for subregion_label in subregion_labels['background_objects']
                            ])
            bootstrap_background_objects=background_objects[bootstrap_sample]
            bootstrap_weights=weights[bootstrap_sample]
            # compute mean transmitted profile for the jackknifed sample
            mean_flux_k, _ , _ = self.mean_profile_around_objects(
                    bootstrap_foreground_objects, bootstrap_background_objects,
                    weights=bootstrap_weights,
                    min_bin=min_bin,max_bin=max_bin,bins=bins,type=type
                    )
            resampled_mean_flux.append(mean_flux_k)
        # calculate bootstrap error
        resampled_mean_flux=np.array(resampled_mean_flux) # convert to array for easy computation
        #bootstrap_variance=np.var(resampled_mean_flux,axis=0,ddof=1)  # = 1/(N-1) * Sum( [f(x)^k - <f(x)^k>]^2 )
        bootstrap_variance=np.nanvar(resampled_mean_flux,axis=0,ddof=1)  # = 1/(N-1) * Sum( [f(x)^k - <f(x)^k>]^2 )
        if return_covariance_matrix:
            stats=np.ma.array(resampled_mean_flux,mask=np.isnan(resampled_mean_flux))
            stats_mean=np.ma.mean(stats,axis=0)
            mat=stats-stats_mean
            Cov=np.ma.dot( (stats-stats_mean).T, stats-stats_mean )  # = sum (f(x_i)^k-<f(x_i)^k>)*(f(x_j)^k-<f(x_j)^k>) from k=1 ... N
            N=(stats-stats_mean).count(axis=0)
            bootstrap_covariance=Cov.data/(N-1)
            return (bootstrap_variance, theta_bins, bootstrap_covariance)
        else:
            return (bootstrap_variance, theta_bins)

    # Landy-Szalay estimator
    # def Landy_Szalay_estimator(self, data_objects, random_objects,
    #                            min_bin=1*u.arcmin, max_bin=30*u.arcmin,
    #                            bins=10, type='linear-bin'):
        
    #     import sys

    #     def logging(msg):
    #         sys.stdout.write(f"\r{msg}")
    #         sys.stdout.flush()

    #     # angular bins
    #     theta_bins, theta_edges=self.angular_bins(min_bin,max_bin,bins,type=type)
    #     # define data-data, data-random, random-random pair counts
    #     DD=np.zeros(bins)
    #     DR=np.zeros(bins)
    #     RR=np.zeros(bins)
    #     N_data=data_objects.size
    #     N_random=random_objects.size
    #     logging('Computing data-data pair counts ...')
    #     for i in range(N_data):
    #         theta=data_objects[i].separation(data_objects[i+1:])
    #         for n in range(bins):
    #             counts=( (theta>=theta_edges[n]) & (theta< theta_edges[n+1]) )
    #             DD[n]+=np.sum(counts)
    #     logging('Computing data-random pair counts ...')
    #     for i in range(N_data):
    #         theta=data_objects[i].separation(random_objects)
    #         for n in range(bins):
    #             counts=( (theta>=theta_edges[n]) & (theta< theta_edges[n+1]) )
    #             DR[n]+=np.sum(counts)
    #     DR=((N_data)/(N_random))*DR # correction for different number of random vs data
    #     logging('Computing random-random pair counts ...')
    #     for i in range(N_random):
    #         theta=random_objects[i].separation(random_objects[i+1:])
    #         for n in range(bins):
    #             counts=( (theta>=theta_edges[n]) & (theta< theta_edges[n+1]) )
    #             RR[n]+=np.sum(counts)
    #     RR=(N_data*(N_data-1))/(N_random*(N_random-1))*RR # correction for different number of random vs data
    #     # Landy-Szalay estimator
    #     #ACF=(DD-2*DR+RR)/RR
    #     ACF = DD/RR-1
    #     ERROR = 1/np.sqrt(DD) 
    #     # integral constraint

    #     return ( ACF, ERROR, theta_bins )
    

    def Landy_Szalay_estimator(self, data_objects, random_objects,
                           min_bin=1*u.arcmin, max_bin=30*u.arcmin,
                           bins=10, type='linear-bin'):

        # angular bins
        theta_bins, theta_edges = self.angular_bins(min_bin, max_bin, bins, type=type)
        N_data = len(data_objects)
        N_random = len(random_objects)

        # --- DD counts ---
        idx1, idx2, sep, _ = data_objects.search_around_sky(data_objects, max_bin)
        # remove self-matches and double counting (i<j)
        mask = idx1 < idx2
        sep = sep[mask]
        DD, _ = np.histogram(sep, bins=theta_edges)

        # --- DR counts ---
        idx1, idx2, sep, _ = data_objects.search_around_sky(random_objects, max_bin)
        DR, _ = np.histogram(sep, bins=theta_edges)
        DR = (N_data / N_random) * DR

        # --- RR counts ---
        idx1, idx2, sep, _ = random_objects.search_around_sky(random_objects, max_bin)
        mask = idx1 < idx2
        sep = sep[mask]
        RR, _ = np.histogram(sep, bins=theta_edges)
        RR = (N_data * (N_data - 1)) / (N_random * (N_random - 1)) * RR

        # --- Landy-Szalay estimator ---
        #ACF = (DD - 2*DR + RR) / RR
        ACF = DD / RR - 1
        ERROR = 1 / np.sqrt(DD)  #  DD.clip(min=1) avoid div-by-zero

        return ACF, ERROR, theta_bins

    # Landy-Szalay estimator: cross
    # def Landy_Szalay_estimator_cross(self,
    #                            data_objects_1, data_objects_2,
    #                            random_objects_1, random_objects_2,
    #                            min_bin=1*u.arcmin, max_bin=30*u.arcmin,
    #                            bins=10, type='linear-bin'):
    #     # angular bins
    #     theta_bins, theta_edges=self.angular_bins(min_bin,max_bin,bins,type=type)
    #     # define data-data, data-random, random-random pair counts
    #     D1D2=np.zeros(bins)
    #     D1R2=np.zeros(bins)
    #     D2R1=np.zeros(bins)
    #     R1R2=np.zeros(bins)
    #     N_data_1=data_objects_1.size
    #     N_data_2=data_objects_2.size
    #     N_random_1=random_objects_1.size
    #     N_random_2=random_objects_2.size
    #     print('Computing data 1 -data 2 pair counts ...')
    #     for i in range(N_data_1):
    #         theta=data_objects_1[i].separation(data_objects_2)
    #         for n in range(bins):
    #             counts=( (theta>=theta_edges[n]) & (theta< theta_edges[n+1]) )
    #             D1D2[n]+=np.sum(counts)
    #     D1D2=D1D2/(N_data_1*N_data_2)
    #     print('Computing data 1 - random 2 pair counts ...')
    #     for i in range(N_data_1):
    #         theta=data_objects_1[i].separation(random_objects_2)
    #         for n in range(bins):
    #             counts=( (theta>=theta_edges[n]) & (theta< theta_edges[n+1]) )
    #             D1R2[n]+=np.sum(counts)
    #     D1R2=D1R2/(N_data_1*N_random_2) # correction for different number of random vs data
    #     print('Computing data 2 - random 1 pair counts ...')
    #     for i in range(N_data_2):
    #         theta=data_objects_2[i].separation(random_objects_1)
    #         for n in range(bins):
    #             counts=( (theta>=theta_edges[n]) & (theta< theta_edges[n+1]) )
    #             D2R1[n]+=np.sum(counts)
    #     D2R1=D2R1/(N_data_2*N_random_1) # correction for different number of random vs data
    #     print('Computing random 1 -random 2 pair counts ...')
    #     for i in range(N_random_1):
    #         theta=random_objects_1[i].separation(random_objects_2)
    #         for n in range(bins):
    #             counts=( (theta>=theta_edges[n]) & (theta< theta_edges[n+1]) )
    #             R1R2[n]+=np.sum(counts)
    #     R1R2=R1R2/(N_random_1*N_random_2) # correction for different number of random vs data
    #     # Landy-Szalay estimator
    #     CCF=(D1D2-D1R2-D2R1+R1R2)/R1R2
    #     ERROR = 1/np.sqrt(D1D2*(N_data_1*N_data_2))
    #     return ( CCF, ERROR, theta_bins )
    
    

    def Landy_Szalay_estimator_cross(self,
                                    data_objects_1, data_objects_2,
                                    random_objects_1, random_objects_2,
                                    min_bin=1*u.arcmin, max_bin=30*u.arcmin,
                                    bins=10, type='linear-bin'):

        # angular bins
        theta_bins, theta_edges = self.angular_bins(min_bin, max_bin, bins, type=type)

        # number counts
        N_data_1 = len(data_objects_1)
        N_data_2 = len(data_objects_2)
        N_random_1 = len(random_objects_1)
        N_random_2 = len(random_objects_2)

        import sys

        def logging(msg):
            sys.stdout.write(f"\r{msg}")
            sys.stdout.flush()

        def pair_counts(cat1, cat2, norm):
            """Compute normalized pair counts using KD-tree search."""
            idx1, idx2, sep, _ = cat1.search_around_sky(cat2, max_bin)
            # sep is already an Angle array
            hist, _ = np.histogram(sep, bins=theta_edges)
            return hist / norm

        logging("Computing data1-data2 pair counts ...")
        D1D2 = pair_counts(data_objects_1, data_objects_2, N_data_1 * N_data_2)

        logging("Computing data1-random2 pair counts ...")
        D1R2 = pair_counts(data_objects_1, random_objects_2, N_data_1 * N_random_2)

        logging("Computing data2-random1 pair counts ...")
        D2R1 = pair_counts(data_objects_2, random_objects_1, N_data_2 * N_random_1)

        logging("Computing random1-random2 pair counts ...")
        R1R2 = pair_counts(random_objects_1, random_objects_2, N_random_1 * N_random_2)

        # Landy-Szalay estimator
        CCF = (D1D2 - D1R2 - D2R1 + R1R2) / R1R2

        # Simple Poisson error estimate
        ERROR = 1 / np.sqrt(D1D2 * N_data_1 * N_data_2)

        return CCF, ERROR, theta_bins
    




    # Davis-Peebles estimator DD/DR-1: cross
    # def Davis_Peebles_estimator_cross(self,
    #                                   data_objects_1, data_objects_2,
    #                                   random_objects_1,
    #                                   min_bin=1*u.arcmin, max_bin=30*u.arcmin,
    #                                   bins=10, type='linear-bin'):
    #     # angular bins
    #     theta_bins, theta_edges=self.angular_bins(min_bin,max_bin,bins,type=type)
    #     # define data-data, data-random, random-random pair counts
    #     D1D2=np.zeros(bins)
    #     D2R1=np.zeros(bins)
    #     N_data_1=data_objects_1.size
    #     N_data_2=data_objects_2.size
    #     N_random_1=random_objects_1.size
    #     print('Computing data 1 - data 2 pair counts ...')
    #     for i in range(N_data_1):
    #         theta=data_objects_1[i].separation(data_objects_2)
    #         for n in range(bins):
    #             counts=( (theta>=theta_edges[n]) & (theta< theta_edges[n+1]) )
    #             D1D2[n]+=np.sum(counts)
    #     D1D2=D1D2/(N_data_1*N_data_2)
    #     print('Computing data 2 - random 1 pair counts ...')
    #     for i in range(N_data_2):
    #         theta=data_objects_2[i].separation(random_objects_1)
    #         for n in range(bins):
    #             counts=( (theta>=theta_edges[n]) & (theta< theta_edges[n+1]) )
    #             D2R1[n]+=np.sum(counts)
    #     D2R1=D2R1/(N_data_2*N_random_1) # correction for different number of random vs data
    #     # Landy-Szalay estimator
    #     CCF=D1D2/D2R1-1

    #     return ( CCF, theta_bins )
    

    def Davis_Peebles_estimator_cross(self,
                                    data_objects_1, data_objects_2,
                                    random_objects_1,
                                    min_bin=1*u.arcmin, max_bin=30*u.arcmin,
                                    bins=10, type='linear-bin'):
        """
        Davis-Peebles estimator for cross-correlation:
                w(theta) = DD / DR - 1
        """

        # angular bins
        theta_bins, theta_edges = self.angular_bins(min_bin, max_bin, bins, type=type)

        N_data_1 = len(data_objects_1)
        N_data_2 = len(data_objects_2)
        N_random_1 = len(random_objects_1)
        
        import sys

        def logging(msg):
            sys.stdout.write(f"\r{msg}")
            sys.stdout.flush()

        def pair_counts(cat1, cat2, norm):
            """KD-tree accelerated pair counting."""
            idx1, idx2, sep, _ = cat1.search_around_sky(cat2, max_bin)
            hist, _ = np.histogram(sep, bins=theta_edges)
            return hist / norm

        logging("Computing data1 – data2 pair counts ...")
        D1D2 = pair_counts(data_objects_1, data_objects_2, N_data_1 * N_data_2)

        logging("Computing data2 – random1 pair counts ...")
        D2R1 = pair_counts(data_objects_2, random_objects_1, N_data_2 * N_random_1)

        # Davis-Peebles estimator
        CCF = D1D2 / D2R1 - 1

        ERROR = 1 / np.sqrt(D1D2 * N_data_1 * N_data_2)

        return CCF, ERROR, theta_bins


    def covariance_matrix(self,stats,type=None):
        """
        Retrun the covarinace matrix of summary statisitcs
            Input:
            stats  : numpy.array, array of resampled statistics
            type   : 'Jackknife' or 'Bootstrap'
        """
        stats=np.ma.array(stats,mask=np.isnan(stats))
        stats_mean=np.ma.mean(stats,axis=0)
        mat=stats-stats_mean
        # Cov = sum (f(x_i)^k-<f(x_i)^k>)*(f(x_j)^k-<f(x_j)^k>) from k=1 ... N
        Cov=np.ma.dot( (stats-stats_mean).T, stats-stats_mean )
        N=(stats-stats_mean).count(axis=0)
        if type=='default':
            return Cov.data
        if type=='Bootstrap':
            return 1/(N-1)*Cov.data
        if type=='Jackknife':
            return (N-1)/N*Cov.data

    def jackknife_error_estimate_Landy_Szalay(self,
                    foreground_objects, background_objects, weights=None,
                    random_objects=None, internal_random=True,
                    min_bin=1*u.arcmin,max_bin=30*u.arcmin,bins=10,
                    jk_region_size=20,plot_check=True,type='linear-bin',
                    return_covariance_matrix=False):
        print('Jackknife error estimation (cross-correlation function)...')
        # initialize random seed
        np.random.seed(seed=0)
        # define jackknife regions from k-mean clustering of random objects
        from kmeans_radec import KMeans, kmeans_sample
        if internal_random: random_objects=self.make_random_catalog(size=50000)
        R=np.vstack( [random_objects.ra.value, random_objects.dec.value] ).T
        km=kmeans_sample(R, jk_region_size, maxiter=100, tol=1.0e-5,verbose=0)
        if plot_check:
            fig,ax=plt.subplots(figsize=(6,5))
            labels=km.find_nearest(R)
            ax.scatter(R[:,0],R[:,1],c=labels,s=1,marker='.')
            ax.set_aspect('equal')
            ax.set_xlabel('RA [deg]')
            ax.set_ylabel('DEC [deg]')
            plt.tight_layout()
            plt.pause(0.05) #plt.show()
        # get the Jackknife labels for all objects
        jk_regions={}
        jk_regions['foreground_objects']=km.find_nearest(
                                            np.vstack( [ foreground_objects.ra.value,
                                                         foreground_objects.dec.value ] ).T
                                                        )
        jk_regions['background_objects']=km.find_nearest(
                                            np.vstack( [ background_objects.ra.value,
                                                         background_objects.dec.value ] ).T
                                                        )
        # angular bins
        theta_bins, theta_edges=self.angular_bins(min_bin,max_bin,bins,type=type)
        # define mean transmitted flux profiles for jackknifed sample
        jk_mean_flux=[]
        for k in range(jk_region_size):
            print('At ',k,'th Jackknife sample ...')
            # jackknife the foreground objects
            jackknife_sample = (jk_regions['foreground_objects'] != k)
            jk_foreground_objects=foreground_objects[jackknife_sample]
            # jackknife the background objects
            jackknife_sample = (jk_regions['background_objects'] != k)
            jk_background_objects=background_objects[jackknife_sample]
            # compute mean transmitted profile for the jackknifed sample
            mean_flux_k, _ = self.Landy_Szalay_estimator(
                    jk_foreground_objects, jk_background_objects,
                    min_bin=min_bin,max_bin=max_bin,bins=bins,type=type
                    )
            jk_mean_flux.append(mean_flux_k)
        # calculate jackknife error
        # np.var(x) is the average of the squared deviations from the mean,
        #   i.e., var = mean(x), where x = abs(a - a.mean())**2.
        jk_mean_flux=np.array(jk_mean_flux) # convert to array for easy computation
        var=np.var(jk_mean_flux,axis=0)  # = (1/N_JK) * Sum( [f(x)^k - <f(x)^k>]^2 )
        variance=(jk_region_size-1)*var # = ((N_JK-1)/N_JK) * Sum( [f(x)^k - <f(x)^k>]^2 )

        if return_covariance_matrix:
            covariance=self.covariance_matrix(jk_mean_flux,type='Jackknife')
            return (variance, theta_bins, covariance)
        else:
            return (variance, theta_bins)

    def jackknife_error_estimate_Landy_Szalay_cross(self,
                    data_objects_1, data_objects_2,
                    random_objects_1, random_objects_2, weights=None,
                    internal_random_objects=None, internal_random=True,
                    min_bin=1*u.arcmin,max_bin=30*u.arcmin,bins=10,
                    jk_region_size=20,plot_check=True,type='linear-bin',
                    return_covariance_matrix=False):
        print('Jackknife error estimation (Landy-Szalay cross-correlation function)...')
        # initialize random seed
        np.random.seed(seed=0)
        # define jackknife regions from k-mean clustering of random objects
        from kmeans_radec import KMeans, kmeans_sample
        if internal_random: internal_random_objects=self.make_random_catalog(size=50000)
        R=np.vstack( [internal_random_objects.ra.value, internal_random_objects.dec.value] ).T
        km=kmeans_sample(R, jk_region_size, maxiter=100, tol=1.0e-5,verbose=0)
        if plot_check:
            fig,ax=plt.subplots(figsize=(6,5))
            labels=km.find_nearest(R)
            ax.scatter(R[:,0],R[:,1],c=labels,s=1,marker='.')
            ax.set_aspect('equal')
            ax.set_xlabel('RA [deg]')
            ax.set_ylabel('DEC [deg]')
            plt.tight_layout()
            plt.pause(0.05) #plt.show()
        # get the Jackknife labels for all objects
        jk_regions={}
        jk_regions['data_objects_1']=km.find_nearest(
                                        np.vstack( [ data_objects_1.ra.value,
                                                     data_objects_1.dec.value ] ).T
                                                    )
        jk_regions['data_objects_2']=km.find_nearest(
                                        np.vstack( [ data_objects_2.ra.value,
                                                     data_objects_2.dec.value ] ).T
                                                    )
        jk_regions['random_objects_1']=km.find_nearest(
                                        np.vstack( [ random_objects_1.ra.value,
                                                     random_objects_1.dec.value ] ).T
                                                    )
        jk_regions['random_objects_2']=km.find_nearest(
                                        np.vstack( [ random_objects_2.ra.value,
                                                     random_objects_2.dec.value ] ).T
                                                    )
        # angular bins
        theta_bins, theta_edges=self.angular_bins(min_bin,max_bin,bins,type=type)
        # define mean transmitted flux profiles for jackknifed sample
        jk_mean_flux=[]
        for k in range(jk_region_size):
            print('At ',k,'th Jackknife sample ...')
            # jackknife the data objects 1 & 2
            jackknife_sample = (jk_regions['data_objects_1'] != k)
            jk_data_objects_1=data_objects_1[jackknife_sample]
            jackknife_sample = (jk_regions['data_objects_2'] != k)
            jk_data_objects_2=data_objects_2[jackknife_sample]
            # jackknife the random objects
            jackknife_sample = (jk_regions['random_objects_1'] != k)
            jk_random_objects_1=random_objects_1[jackknife_sample]
            jackknife_sample = (jk_regions['random_objects_2'] != k)
            jk_random_objects_2=random_objects_2[jackknife_sample]
            # compute mean transmitted profile for the jackknifed sample
            mean_flux_k, _ = self.Landy_Szalay_estimator_cross(
                    jk_data_objects_1, jk_data_objects_2,
                    jk_random_objects_1, jk_random_objects_2,
                    min_bin=min_bin,max_bin=max_bin,bins=bins,type=type
                    )
            jk_mean_flux.append(mean_flux_k)
        # calculate jackknife error
        # np.var(x) is the average of the squared deviations from the mean,
        #   i.e., var = mean(x), where x = abs(a - a.mean())**2.
        jk_mean_flux=np.array(jk_mean_flux) # convert to array for easy computation
        var=np.var(jk_mean_flux,axis=0)  # = (1/N_JK) * Sum( [f(x)^k - <f(x)^k>]^2 )
        variance=(jk_region_size-1)*var # = ((N_JK-1)/N_JK) * Sum( [f(x)^k - <f(x)^k>]^2 )

        if return_covariance_matrix:
            covariance=self.covariance_matrix(jk_mean_flux,type='Jackknife')
            return (variance, theta_bins, covariance)
        else:
            return (variance, theta_bins)

    def jackknife_error_estimate_Davis_Peebles_cross(self,
                    data_objects_1, data_objects_2,
                    random_objects_1, weights=None,
                    internal_random_objects=None, internal_random=True,
                    min_bin=1*u.arcmin,max_bin=30*u.arcmin,bins=10,
                    jk_region_size=20,plot_check=True,type='linear-bin',
                    return_covariance_matrix=False):
        print('Jackknife error estimation (Davis-Peebles cross-correlation function)...')
        # initialize random seed
        np.random.seed(seed=0)
        # define jackknife regions from k-mean clustering of random objects
        from kmeans_radec import KMeans, kmeans_sample
        if internal_random: internal_random_objects=self.make_random_catalog(size=50000)
        R=np.vstack( [internal_random_objects.ra.value, internal_random_objects.dec.value] ).T
        km=kmeans_sample(R, jk_region_size, maxiter=100, tol=1.0e-5,verbose=0)
        if plot_check:
            fig,ax=plt.subplots(figsize=(6,5))
            labels=km.find_nearest(R)
            ax.scatter(R[:,0],R[:,1],c=labels,s=1,marker='.')
            ax.set_aspect('equal')
            ax.set_xlabel('RA [deg]')
            ax.set_ylabel('DEC [deg]')
            plt.tight_layout()
            plt.pause(0.05) #plt.show()
        # get the Jackknife labels for all objects
        jk_regions={}
        jk_regions['data_objects_1']=km.find_nearest(
                                        np.vstack( [ data_objects_1.ra.value,
                                                     data_objects_1.dec.value ] ).T
                                                    )
        jk_regions['data_objects_2']=km.find_nearest(
                                        np.vstack( [ data_objects_2.ra.value,
                                                     data_objects_2.dec.value ] ).T
                                                    )
        jk_regions['random_objects_1']=km.find_nearest(
                                        np.vstack( [ random_objects_1.ra.value,
                                                     random_objects_1.dec.value ] ).T
                                                    )
        # angular bins
        theta_bins, theta_edges=self.angular_bins(min_bin,max_bin,bins,type=type)
        # define mean transmitted flux profiles for jackknifed sample
        jk_mean_flux=[]
        for k in range(jk_region_size):
            print('At ',k,'th Jackknife sample ...')
            # jackknife the data objects 1 & 2
            jackknife_sample = (jk_regions['data_objects_1'] != k)
            jk_data_objects_1=data_objects_1[jackknife_sample]
            jackknife_sample = (jk_regions['data_objects_2'] != k)
            jk_data_objects_2=data_objects_2[jackknife_sample]
            # jackknife the random objects
            jackknife_sample = (jk_regions['random_objects_1'] != k)
            jk_random_objects_1=random_objects_1[jackknife_sample]
            # compute mean transmitted profile for the jackknifed sample
            mean_flux_k, _ = self.Davis_Peebles_estimator_cross(
                    jk_data_objects_1, jk_data_objects_2,
                    jk_random_objects_1,
                    min_bin=min_bin,max_bin=max_bin,bins=bins,type=type
                    )
            jk_mean_flux.append(mean_flux_k)
        # calculate jackknife error
        # np.var(x) is the average of the squared deviations from the mean,
        #   i.e., var = mean(x), where x = abs(a - a.mean())**2.
        jk_mean_flux=np.array(jk_mean_flux) # convert to array for easy computation
        var=np.var(jk_mean_flux,axis=0)  # = (1/N_JK) * Sum( [f(x)^k - <f(x)^k>]^2 )
        variance=(jk_region_size-1)*var # = ((N_JK-1)/N_JK) * Sum( [f(x)^k - <f(x)^k>]^2 )

        if return_covariance_matrix:
            covariance=self.covariance_matrix(jk_mean_flux,type='Jackknife')
            return (variance, theta_bins, covariance)
        else:
            return (variance, theta_bins)


    # estimator: Lya forest auto-correlation function
    def auto_correlation_function(self, objects, weights=None,
                                  min_bin=1*u.arcmin,max_bin=30*u.arcmin,bins=10,
                                  type='linear-bin'):
        # angular bins
        theta_bins, theta_edges=self.angular_bins(min_bin,max_bin,bins,type=type)
        xi=np.zeros(bins)
        npairs=np.zeros(bins)
        # copmute auto-correlation function
        print('computing auto-correlation function ...')
        for i in range(0,len(objects)):
            for j in range(i,len(objects)):
                theta=objects[i].separation(objects[j])
                for n in range(bins):
                    if ( (theta>=theta_edges[n]) &
                         (theta< theta_edges[n+1]) ):
                        xi[n]+=weights[i]*weights[j]
                        npairs[n]+=1
        xi=xi/npairs
        return(xi, npairs, theta_bins)

    # estimator: Lya forest auto-correlation function
    def cross_correlation_function(self, objects1, objects2, weights1=None, weights2=None,
                                  min_bin=1*u.arcmin,max_bin=30*u.arcmin,bins=10,
                                  type='linear-bin'):
        # angular bins
        theta_bins, theta_edges=self.angular_bins(min_bin,max_bin,bins,type=type)
        xi=np.zeros(bins)
        npairs=np.zeros(bins)
        # copmute auto-correlation function
        print('computing cross-correlation function ...')
        for i in range(len(objects1)):
            for j in range(len(objects2)):
                theta=objects1[i].separation(objects2[j])
                for n in range(bins):
                    if ( (theta>=theta_edges[n]) &
                         (theta< theta_edges[n+1]) ):
                        xi[n]+=weights1[i]*weights2[j]
                        npairs[n]+=1
        xi=xi/npairs
        return(xi, npairs, theta_bins)
    # Jackknife error estimator: auto-correlation function
    def jackknife_error_estimate_auto(self,objects, weights=None,
                                      random_objects=None, internal_random=True,
                                      min_bin=1*u.arcmin,max_bin=30*u.arcmin,bins=10,
                                      jk_region_size=20,plot_check=True,type='linear-bin',
                                      seed=123456):
        print('Jackknife error estimation (auto-correlation function)...')
        # initialize random seed
        np.random.seed(seed=seed)
        # define jackknife regions from k-mean clustering of random objects
        from kmeans_radec import KMeans, kmeans_sample
        if internal_random: random_objects=self.make_random_catalog(size=50000)
        R=np.vstack( [random_objects.ra.value, random_objects.dec.value] ).T
        km=kmeans_sample(R, jk_region_size, maxiter=100, tol=1.0e-5)
        if plot_check:
            fig,ax=plt.subplots(figsize=(6,5))
            labels=km.find_nearest(R)
            ax.scatter(R[:,0],R[:,1],c=labels,s=1,marker='.')
            ax.set_aspect('equal')
            ax.set_xlabel('RA [deg]')
            ax.set_ylabel('DEC [deg]')
            plt.tight_layout()
            plt.pause(0.05) #plt.show()
        # get the Jackknife labels for all objects
        jk_regions={}
        jk_regions['objects']=km.find_nearest(
                                            np.vstack( [ objects.ra.value,
                                                         objects.dec.value ] ).T
                                             )
        # angular bins
        theta_bins, theta_edges=self.angular_bins(min_bin,max_bin,bins,type=type)
        jk_xi=[]
        for k in range(jk_region_size):
            print('At ',k,'th Jackknife sample ...')
            jackknife_sample = (jk_regions['objects'] != k)
            jk_objects=objects[jackknife_sample]
            jk_weights=weights[jackknife_sample]
            # compute auto-correlation function for jackknifed sample
            xi_k, _ , _ = self.auto_correlation_function(
                                jk_objects, weights=jk_weights,
                                min_bin=min_bin,max_bin=max_bin,bins=bins,
                                type=type
                                )
            jk_xi.append(xi_k)
        # calculate jackknife error
        # np.var(x) is the average of the squared deviations from the mean,
        #   i.e., var = mean(x), where x = abs(a - a.mean())**2.
        jk_xi=np.array(jk_xi) # convert to array for easy computation
        var=np.var(jk_xi,axis=0)  # = (1/N_JK) * Sum( [f(x)^k - <f(x)^k>]^2 )
        jackknife_variance=(jk_region_size-1)*var # = ((N_JK-1)/N_JK) * Sum( [f(x)^k - <f(x)^k>]^2 )
        return (jackknife_variance, theta_bins)
    # Jackknife error estimator: cross-correlation function
    def jackknife_error_estimate_cross(self,objects1, objects2, weights1=None, weights2=None,
                                      random_objects=None, internal_random=True,
                                      min_bin=1*u.arcmin,max_bin=30*u.arcmin,bins=10,
                                      jk_region_size=20,plot_check=True,type='linear-bin',
                                      seed=123456):
        print('Jackknife error estimation (cross-correlation function)...')
        # initialize random seed
        np.random.seed(seed=seed)
        # define jackknife regions from k-mean clustering of random objects
        from kmeans_radec import KMeans, kmeans_sample
        if internal_random: random_objects=self.make_random_catalog(size=50000)
        R=np.vstack( [random_objects.ra.value, random_objects.dec.value] ).T
        km=kmeans_sample(R, jk_region_size, maxiter=100, tol=1.0e-5)
        if plot_check:
            fig,ax=plt.subplots(figsize=(6,5))
            labels=km.find_nearest(R)
            ax.scatter(R[:,0],R[:,1],c=labels,s=1,marker='.')
            ax.set_aspect('equal')
            ax.set_xlabel('RA [deg]')
            ax.set_ylabel('DEC [deg]')
            plt.tight_layout()
            plt.pause(0.05) #plt.show()
        # get the Jackknife labels for all objects
        jk_regions={}
        jk_regions['objects1']=km.find_nearest(
                                            np.vstack( [ objects1.ra.value,
                                                         objects1.dec.value ] ).T
                                             )
        jk_regions['objects2']=km.find_nearest(
                                            np.vstack( [ objects2.ra.value,
                                                         objects2.dec.value ] ).T
                                             )
        # angular bins
        theta_bins, theta_edges=self.angular_bins(min_bin,max_bin,bins,type=type)
        jk_xi=[]
        for k in range(jk_region_size):
            print('At ',k,'th Jackknife sample ...')
            jackknife_sample = (jk_regions['objects1'] != k)
            jk_objects1=objects1[jackknife_sample]
            jk_weights1=weights1[jackknife_sample]

            jackknife_sample = (jk_regions['objects2'] != k)
            jk_objects2=objects2[jackknife_sample]
            jk_weights2=weights2[jackknife_sample]
            # compute cross-correlation function for jackknifed sample
            xi_k, _ , _ = self.cross_correlation_function(
                                jk_objects1, jk_objects2,
                                weights1=jk_weights1, weights2=jk_weights2,
                                min_bin=min_bin,max_bin=max_bin,bins=bins,type=type
                                )
            jk_xi.append(xi_k)
        # calculate jackknife error
        # np.var(x) is the average of the squared deviations from the mean,
        #   i.e., var = mean(x), where x = abs(a - a.mean())**2.
        jk_xi=np.array(jk_xi) # convert to array for easy computation
        var=np.var(jk_xi,axis=0)  # = (1/N_JK) * Sum( [f(x)^k - <f(x)^k>]^2 )
        jackknife_variance=(jk_region_size-1)*var # = ((N_JK-1)/N_JK) * Sum( [f(x)^k - <f(x)^k>]^2 )
        return (jackknife_variance, theta_bins)

    # Jackknife error estimator: mean IGM transmission
    def jackknife_error_estimate_mean(self,
                    background_objects, weights=None,
                    random_objects=None, internal_random=True,
                    jk_region_size=20,plot_check=True):
        print('Jackknife error estimation (mean IGM transmission)...')
        # initialize random seed
        np.random.seed(seed=0)
        # define jackknife regions from k-mean clustering of random objects
        from kmeans_radec import KMeans, kmeans_sample
        if internal_random: random_objects=self.make_random_catalog(size=50000)
        R=np.vstack( [random_objects.ra.value, random_objects.dec.value] ).T
        km=kmeans_sample(R, jk_region_size, maxiter=100, tol=1.0e-5)
        if plot_check:
            fig,ax=plt.subplots(figsize=(6,5))
            labels=km.find_nearest(R)
            ax.scatter(R[:,0],R[:,1],c=labels,s=1,marker='.')
            ax.set_aspect('equal')
            ax.set_xlabel('RA [deg]')
            ax.set_ylabel('DEC [deg]')
            plt.tight_layout()
            plt.pause(0.05) #plt.show()
        # get the Jackknife labels for all objects
        jk_regions={}
        jk_regions['background_objects']=km.find_nearest(
                                            np.vstack( [ background_objects.ra.value,
                                                         background_objects.dec.value ] ).T
                                                        )
        # define mean transmitted flux profiles for jackknifed sample
        jk_mean_flux=[]
        for k in range(jk_region_size):
            print('At ',k,'th Jackknife sample ...')
            # jackknife the background objects
            jackknife_sample_bg = (jk_regions['background_objects'] != k)
            jk_background_objects=background_objects[jackknife_sample_bg]
            jk_weights=weights[jackknife_sample_bg]
            # compute mean IGM transmission for the jackknifed sample
            mean_flux_k=np.mean(jk_weights)
            jk_mean_flux.append(mean_flux_k)
        # calculate jackknife error
        # np.var(x) is the average of the squared deviations from the mean,
        #   i.e., var = mean(x), where x = abs(a - a.mean())**2.
        jk_mean_flux=np.array(jk_mean_flux) # convert to array for easy computation
        var=np.var(jk_mean_flux)  # = (1/N_JK) * Sum( [f(x)^k - <f(x)^k>]^2 )
        jackknife_variance=(jk_region_size-1)*var # = ((N_JK-1)/N_JK) * Sum( [f(x)^k - <f(x)^k>]^2 )
        return (jackknife_variance)

    # Bootstrapping error estimator: mean Lya transmission
    def bootstrap_error_estimate_mean(self,
                    background_objects, weights=None,
                    random_objects=None, internal_random=True,
                    subregion_size=20,resampling_size=20,plot_check=True):
        print('Bootstrapping error estimation (mean IGM transmission)...')
        # initialize random seed
        np.random.seed(seed=0)
        # define subvolume regions from k-mean clustering of random objects
        from kmeans_radec import KMeans, kmeans_sample
        if internal_random: random_objects=self.make_random_catalog(size=50000)
        R=np.vstack( [random_objects.ra.value, random_objects.dec.value] ).T
        km=kmeans_sample(R, subregion_size, maxiter=100, tol=1.0e-5, verbose=0)
        if plot_check:
            fig,ax=plt.subplots(figsize=(6,5))
            labels=km.find_nearest(R)
            ax.scatter(R[:,0],R[:,1],c=labels,s=1,marker='.')
            ax.set_aspect('equal')
            ax.set_xlabel('RA [deg]')
            ax.set_ylabel('DEC [deg]')
            plt.tight_layout()
            plt.pause(0.05) #plt.show()
        # get the subregion labels for all objects
        subregion_labels={}
        subregion_labels['background_objects']=km.find_nearest(
                                                np.vstack( [ background_objects.ra.value,
                                                             background_objects.dec.value ] ).T
                                                             )
        # bootstrap resampling the subregions
        resampled_regions=resampling_size*[None]
        for k in range(resampling_size):
            resampled_regions[k]=np.random.choice(
                                            np.unique(km.labels),
                                            size=subregion_size,
                                            replace=True
                                            )
        # define mean transmitted flux profiles for resampled regions
        resampled_mean_flux=[]
        for k in range(resampling_size):
            print('At ',k,'th Bootstrapping sample ...')
            # booststrap resampling of the background objects
            bootstrap_sample=np.array([
                            (subregion_label in resampled_regions[k])
                            for subregion_label in subregion_labels['background_objects']
                            ])
            bootstrap_background_objects=background_objects[bootstrap_sample]
            bootstrap_weights=weights[bootstrap_sample]
            # compute mean transmitted profile for the jackknifed sample
            mean_flux_k=np.mean(bootstrap_weights)
            resampled_mean_flux.append(mean_flux_k)
        # calculate bootstrap error
        resampled_mean_flux=np.array(resampled_mean_flux) # convert to array for easy computation
        #bootstrap_variance=np.var(resampled_mean_flux,axis=0,ddof=1)  # = 1/(N-1) * Sum( [f(x)^k - <f(x)^k>]^2 )
        bootstrap_variance=np.nanvar(resampled_mean_flux,ddof=1)  # = 1/(N-1) * Sum( [f(x)^k - <f(x)^k>]^2 )
        return bootstrap_variance


# main test
if __name__ == '__main__':
    hsc=HSC_toolkits()
    selected_specz=hsc.select_background_objects(zmin=4.98,zmax=5.89,filter='NB0718',
                                                 specz_flag=True, magcut_flag=False)

    pd.set_option('display.max_rows', hsc.specz.shape[0]+1)
#    pd.set_option('display.max_columns', hsc.specz.shape[1]+1)

    print(selected_specz[ [ 'specz_redshift',
                            'specz_duplicationflag',
                            'specz_flag_3dhst_v4_1_5',
                            #'specz_flag_sdss_dr15',  # All False
                            #'specz_flag_zcosmos_bright_dr3', # All False
                            #'specz_flag_gama_dr3', # All False
                            #'specz_flag_udsz_dr1', # All False
                            #'specz_flag_vandels_dr2', # All False
                            #'specz_flag_c3r2_dr2', # All False
                            #'specz_flag_vvds_drfinal', # All False
                            'specz_flag_deimos_2018', # All False
                            #'specz_flag_fmos_dr2', # All False
                            #'specz_flag_lega_c_dr2', # All False
                            #'specz_flag_primus_dr1', # All False
                            #'specz_flag_vipers_dr2', # All False
                            #'specz_flag_wigglez_dr1', # All False
                            #'specz_flag_deep23_dr4_egs', # All False
                            #'specz_flag_2dfgrs', # All False
                            #'specz_flag_6dfgrs' # All False
                            ] ] )

    print(selected_specz[['ra','dec','tract','patch_s','z_apertureflux_15_flux','z_apertureflux_15_mag','z_apertureflux_15_magerr']])

    print(selected_specz[['ra','dec','specz_redshift']])

    fig,ax = plt.subplots(figsize=(5.3, 4))
    #plt.rc('text', usetex=True)
    plt.rc('font', family='Times New Roman')
    matplotlib.rcParams.update({'font.size': 16})

    ax.hist(selected_specz['z_apertureflux_15_mag'],range=(22,29),bins=14,edgecolor='black', linewidth=1, label='All')

    is_3dhst = selected_specz['specz_flag_3dhst_v4_1_5']==True
    ax.hist(selected_specz['z_apertureflux_15_mag'][is_3dhst],range=(22,29),bins=14,color='green',edgecolor='black', linewidth=1,
            alpha=0.5,label='3D-HST')

    is_deimos10k = selected_specz['specz_flag_deimos_2018']==True
    ax.hist(selected_specz['z_apertureflux_15_mag'][is_deimos10k],range=(22,29),bins=14,color='red',edgecolor='black', linewidth=1,
            alpha=0.5,label='DEIMOS 10k')


    ax.legend(frameon=False)
    ax.set_xlabel('z mag (1.5\") [mag]')
    ax.set_xlim(22,29)
    ax.set_xticks(arange(22,29,1.0))
    ax.set_ylabel('Number of background sources')
    tight_layout()
    show()
