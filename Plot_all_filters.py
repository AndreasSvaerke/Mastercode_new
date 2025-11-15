import numpy as np
import matplotlib.pyplot as plt
from astropy.wcs import WCS
import os

def nan_mean_std(data):
    m = np.nanmean(data)
    s = np.nanstd(data)
    return m,s

def find_coor_in_ICRS(header,x,y):
    w = WCS(header)
    sky = w.pixel_to_world(x,y)   
    return sky

def find_coor_in_physical(header,sky):
    w = WCS(header)
    x,y = w.world_to_pixel(sky)
    return x,y


def find_diff_bands_cutout_save(output_folder, obs, data, header, skyposition, ep = 1, inter = 25):
    '''
    Takes the data and headers from all the filters.
    The skypostions are from the F115W filters difference image
    obs is a string with the ID (fx. 'obs17')
    Saves the cutouts of the images in the output folder.
    '''
    mean, std = nan_mean_std(data[2])

    for i in range(len(skyposition)):
        sky = skyposition[i]

        old_f115w_x, old_f115w_y = find_coor_in_physical(header[0],sky)
        new_f115w_x, new_f115w_y = find_coor_in_physical(header[1],sky)
        diff_f115w_x, diff_f115w_y = find_coor_in_physical(header[2],sky)
        f150w_x, f150w_y = find_coor_in_physical(header[3],sky)
        f200w_x, f200w_y = find_coor_in_physical(header[4],sky)
        f277w_x, f277w_y = find_coor_in_physical(header[5],sky)
        f356w_x, f356w_y = find_coor_in_physical(header[6],sky)
        f444w_x, f444w_y = find_coor_in_physical(header[7],sky)
        
        diff_f115w_x, diff_f115w_y = int(diff_f115w_x), int(diff_f115w_y)
        new_f115w_x, new_f115w_y = int(new_f115w_x), int(new_f115w_y)
        old_f115w_x, old_f115w_y = int(old_f115w_x), int(old_f115w_y)
        f150w_x, f150w_y = int(f150w_x), int(f150w_y)
        f200w_x, f200w_y = int(f200w_x), int(f200w_y)
        f277w_x, f277w_y = int(f277w_x), int(f277w_y)
        f356w_x, f356w_y = int(f356w_x), int(f356w_y)
        f444w_x, f444w_y = int(f444w_x), int(f444w_y)

        fig, axs = plt.subplots(2, 4, figsize=(16, 10))
        # Plot data on each subplot

        axs[0,0].imshow(data[0][old_f115w_y-inter:old_f115w_y+inter, old_f115w_x-inter:old_f115w_x+inter], vmin = mean-ep*std, vmax = mean+ep*std)
        axs[0,0].set_title('Old_F115W')
        #axs[0,0].axis("off")
        
        axs[0,1].imshow(data[1][new_f115w_y-inter:new_f115w_y+inter, new_f115w_x-inter:new_f115w_x+inter], vmin = mean-ep*std, vmax = mean+ep*std)
        axs[0,1].set_title('New_F115W')
        #axs[0,1].axis("off")
         
        axs[0,2].imshow(data[2][diff_f115w_y-inter:diff_f115w_y+inter, diff_f115w_x-inter:diff_f115w_x+inter], vmin = mean-ep*std, vmax = mean+ep*std)
        axs[0,2].set_title('Diff_F115W')
        #axs[0,2].axis("off")

        axs[0,3].imshow(data[3][f150w_y-inter:f150w_y+inter, f150w_x-inter:f150w_x+inter], vmin = mean-ep*std, vmax = mean+ep*std)
        axs[0,3].set_title('F150W')
        #axs[0,3].axis("off")

        axs[1,0].imshow(data[4][f200w_y-inter:f200w_y+inter, f200w_x-inter:f200w_x+inter], vmin = mean-ep*std, vmax = mean+ep*std)
        axs[1,0].set_title('F200W')
        #axs[1,0].axis("off")
        
        axs[1,1].imshow(data[5][f277w_y-inter:f277w_y+inter, f277w_x-inter:f277w_x+inter], vmin = mean-ep*std, vmax = mean+ep*std)
        axs[1,1].set_title('F277W')
        #axs[1,1].axis("off")
        
        axs[1,2].imshow(data[6][f356w_y-inter:f356w_y+inter, f356w_x-inter:f356w_x+inter], vmin = mean-ep*std, vmax = mean+ep*std)
        axs[1,2].set_title('F356W')
        #axs[1,2].axis("off")

        axs[1,3].imshow(data[7][f444w_y-inter:f444w_y+inter, f444w_x-inter:f444w_x+inter], vmin = mean-ep*std, vmax = mean+ep*std)
        axs[1,3].set_title('F444W')
        #axs[1,3].axis("off")
        
        
        # Adjust layout to prevent overlap
        plt.tight_layout()
       
        ra,dec = sky.ra.value,sky.dec.value

        ra_str = str(ra)
        dec_str = str(dec)

        # Save each figure with a unique name
        filename = os.path.join(output_folder, f"{obs}_{ra_str}_{dec_str}.png")
        plt.savefig(filename, dpi=300)
        plt.close(fig)  # Close the figure to free memory

        # Show the plot
        plt.show()


def find_diff_bands_cutout_plot(obs, data, header, skyposition, ep = 1, inter = 25):
    '''
    Takes the data and headers from all the filters.
    The skypostions are from the F115W filters difference image
    obs is a string with the ID (fx. 'obs17')
    Saves the cutouts of the images in the output folder.
    '''
    mean, std = nan_mean_std(data[2])

    for i in range(len(skyposition)):
        sky = skyposition[i]

        old_f115w_x, old_f115w_y = find_coor_in_physical(header[0],sky)
        new_f115w_x, new_f115w_y = find_coor_in_physical(header[1],sky)
        diff_f115w_x, diff_f115w_y = find_coor_in_physical(header[2],sky)
        f150w_x, f150w_y = find_coor_in_physical(header[3],sky)
        f200w_x, f200w_y = find_coor_in_physical(header[4],sky)
        f277w_x, f277w_y = find_coor_in_physical(header[5],sky)
        f356w_x, f356w_y = find_coor_in_physical(header[6],sky)
        f444w_x, f444w_y = find_coor_in_physical(header[7],sky)
        
        diff_f115w_x, diff_f115w_y = int(diff_f115w_x), int(diff_f115w_y)
        new_f115w_x, new_f115w_y = int(new_f115w_x), int(new_f115w_y)
        old_f115w_x, old_f115w_y = int(old_f115w_x), int(old_f115w_y)
        f150w_x, f150w_y = int(f150w_x), int(f150w_y)
        f200w_x, f200w_y = int(f200w_x), int(f200w_y)
        f277w_x, f277w_y = int(f277w_x), int(f277w_y)
        f356w_x, f356w_y = int(f356w_x), int(f356w_y)
        f444w_x, f444w_y = int(f444w_x), int(f444w_y)

        fig, axs = plt.subplots(2, 4, figsize=(16, 10))
        # Plot data on each subplot

        axs[0,0].imshow(data[0][old_f115w_y-inter:old_f115w_y+inter, old_f115w_x-inter:old_f115w_x+inter], vmin = mean-ep*std, vmax = mean+ep*std)
        axs[0,0].set_title('Old_F115W')
        #axs[0,0].axis("off")
        
        axs[0,1].imshow(data[1][new_f115w_y-inter:new_f115w_y+inter, new_f115w_x-inter:new_f115w_x+inter], vmin = mean-ep*std, vmax = mean+ep*std)
        axs[0,1].set_title('New_F115W')
        #axs[0,1].axis("off")
         
        axs[0,2].imshow(data[2][diff_f115w_y-inter:diff_f115w_y+inter, diff_f115w_x-inter:diff_f115w_x+inter], vmin = mean-ep*std, vmax = mean+ep*std)
        axs[0,2].set_title('Diff_F115W')
        #axs[0,2].axis("off")

        axs[0,3].imshow(data[3][f150w_y-inter:f150w_y+inter, f150w_x-inter:f150w_x+inter], vmin = mean-ep*std, vmax = mean+ep*std)
        axs[0,3].set_title('F150W')
        #axs[0,3].axis("off")

        axs[1,0].imshow(data[4][f200w_y-inter:f200w_y+inter, f200w_x-inter:f200w_x+inter], vmin = mean-ep*std, vmax = mean+ep*std)
        axs[1,0].set_title('F200W')
        #axs[1,0].axis("off")
        
        axs[1,1].imshow(data[5][f277w_y-inter:f277w_y+inter, f277w_x-inter:f277w_x+inter], vmin = mean-ep*std, vmax = mean+ep*std)
        axs[1,1].set_title('F277W')
        #axs[1,1].axis("off")
        
        axs[1,2].imshow(data[6][f356w_y-inter:f356w_y+inter, f356w_x-inter:f356w_x+inter], vmin = mean-ep*std, vmax = mean+ep*std)
        axs[1,2].set_title('F356W')
        #axs[1,2].axis("off")

        axs[1,3].imshow(data[7][f444w_y-inter:f444w_y+inter, f444w_x-inter:f444w_x+inter], vmin = mean-ep*std, vmax = mean+ep*std)
        axs[1,3].set_title('F444W')
        #axs[1,3].axis("off")
        
        
        # Adjust layout to prevent overlap
        plt.tight_layout()
       
        ra,dec = sky.ra.value,sky.dec.value

        ra_str = str(ra)
        dec_str = str(dec)

        print(f"{obs} RA: {ra_str}, DEC: {dec_str} i: {i}")

        # Show the plot
        plt.show()