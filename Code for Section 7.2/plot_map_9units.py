import numpy as np
import pandas as pd
import shapefile as shp
import matplotlib.pyplot as plt
import matplotlib.colorbar as cbar
from colour import Color
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import seaborn as sns
import pickle
import os


# Converting Shapefile Data Into Pandas Dataframes:
def read_shapefile(sf):
    #fetching the headings from the shape file
    fields = [x[0] for x in sf.fields][1:]
    #fetching the records from the shape file
    records = [list(i) for i in sf.records()]
    shps = [s.points for s in sf.shapes()]
    #converting shapefile data into pandas dataframe
    df = pd.DataFrame(columns=fields, data=records)
    #assigning the coordinates
    df = df.assign(coords=shps)
    return df


def plot_map(df, text = True, x_lim = None, y_lim = None, id_list = None, figsize = (11,9)):
    plt.figure(figsize = figsize)
    id=0
    for coord in df.coords:
        # plot points
        x = [i[0] for i in coord]
        y = [i[1] for i in coord]
        plt.plot(x, y, 'k')
        # 
        if text & (x_lim == None) & (y_lim == None):
            x0 = np.mean(x)
            y0 = np.mean(y)
            plt.text(x0, y0, id, fontsize=10)
        id = id+1
    
    if (x_lim != None) & (y_lim != None):     
        plt.xlim(x_lim)
        plt.ylim(y_lim)

def calc_color(data, color, num_color, bins):
    # unique_values = len(set(data))
    # num_color = min(num_color, unique_values - 1)
    if color   == 1:
        color_sq =  ['#dadaebFF','#bcbddcF0','#9e9ac8F0','#807dbaF0','#6a51a3F0','#54278fF0'];
        color_sq = list(Color("#dadaeb").range_to(Color("#6a51a3"),num_color))
        color_sq = [c.hex for c in color_sq]
        colors = 'Purples';
    elif color == 2:
        color_sq = ['#effcef','#94d3ac','#50d890','#55ae95'];
        colors = 'Greens';
        color_sq = list(Color(color_sq[0]).range_to(Color(color_sq[-1]),num_color))
        color_sq = [c.hex for c in color_sq]
    elif color == 3:
        color_sq = ['#f7f7f7','#d9d9d9','#bdbdbd','#969696','#636363','#252525'];
        colors = 'Greys';
        color_sq = list(Color(color_sq[0]).range_to(Color(color_sq[-1]),num_color))
        color_sq = [c.hex for c in color_sq]
    elif color == 4:
        color_sq = ['#f8fbff','#eaf4ff', '#d6eaff','#add6ff','#84c1ff'];
        colors = 'Blues';
        color_sq = list(Color(color_sq[0]).range_to(Color(color_sq[-1]),num_color))
        color_sq = [c.hex for c in color_sq]
    elif color == 5:  # Pink color option
        color_sq = ['#ffe6f0', '#ffb3d9', '#ff80bf', '#ff4da6', '#ff1a8c']
        color_sq = list(Color(color_sq[0]).range_to(Color(color_sq[-1]), num_color))
        color_sq = [c.hex for c in color_sq]
        colors = color_sq  # Directly use the color
    elif color == 6:
        color_sq = ['#FDF8E1', '#FCEFB4', '#FAE588', '#F9DC5C', '#F9DC5C']
        color_sq = list(Color(color_sq[0]).range_to(Color(color_sq[-1]), num_color))
        color_sq = [c.hex for c in color_sq]
        colors = color_sq  # Directly use the color
    elif color == 7:
        color_sq = ['#FFE5EC', '#FFC2D1', '#FFB3C6', '#FF8FAB', '#FB6F92']
        color_sq = list(Color(color_sq[0]).range_to(Color(color_sq[-1]), num_color))
        color_sq = [c.hex for c in color_sq]
        colors = color_sq  # Directly use the color
    elif color == 8:
        color_sq = ['#C7F9CC', '#80ED99', '#57CC99', '#38A3A5', '#22577A']
        color_sq = list(Color(color_sq[0]).range_to(Color(color_sq[-1]), num_color))
        color_sq = [c.hex for c in color_sq]
        colors = color_sq  # Directly use the color
    elif color == 9:
        color_sq = ['#ff0000','#ff0000','#ff0000','#ff0000','#ff0000','#ff0000'];
        color_sq = list(Color(color_sq[0]).range_to(Color(color_sq[-1]),num_color))
        color_sq = [c.hex for c in color_sq]
    else:
        color_sq = ['#ffffd4','#fee391','#fec44f','#fe9929','#d95f0e','#993404'];
        colors = 'YlOrBr';
        num_color = 15
        color_sq = list(Color(color_sq[0]).range_to(Color(color_sq[-1]),num_color))
        color_sq = [c.hex for c in color_sq]

    if len(bins) == 0:
        new_data, bins = pd.qcut(data, num_color, retbins=True, labels=False, duplicates='drop')
    else:
        new_data, bins = pd.cut(data, bins=bins, labels=False, include_lowest=True, retbins=True)
    if len(bins) - 1 != len(color_sq):
        color_sq = list(Color(color_sq[0]).range_to(Color(color_sq[-1]), len(bins) - 1))
        color_sq = [c.hex for c in color_sq]
    color_ton = [color_sq[val] for val in new_data]
    # color_ton = []
    # for val in new_data:
    #     color_ton.append(color_sq[val])
    if color != 9:
        colors = sns.color_palette(colors, n_colors=num_color)
        sns.palplot(colors, 0.6);
    return color_ton, bins, color_sq

def plot_cities_data(df, title, ids, data=None,color=None, num_color=15, print_id=False, bins=[]):
    color_ton, bins, color_sq = calc_color(data, color, num_color, bins)
    ax = plot_map_fill_multiples_ids_tone(df, title, ids, print_id, color_ton, bins, data, color_sq, x_lim = None, y_lim = None);
    return ax

def plot_map_fill_multiples_ids_tone(df, title, ids,  
                                     print_id, color_ton, 
                                     bins, data = None,
                                     color_sq = None,
                                     x_lim = None, 
                                     y_lim = None, 
                                     show_title = False,
                                     figsize = (12,9)):
   
        
    plt.figure(figsize = figsize)
    fig, ax = plt.subplots(figsize = figsize)
    if show_title:
        fig.suptitle(title, fontsize=16)
    for coord in df.coords:
        x = [i[0] for i in coord]
        y = [i[1] for i in coord]
        ax.plot(x, y, 'k')
            
    for id in ids:
        #shape_ex = sf.shape(id)
        points = df.coords[id]
        x_lon = np.zeros((len(points),1))
        y_lat = np.zeros((len(points),1))
        for ip in range(len(points)):
            x_lon[ip], y_lat[ip] = points[ip]
        ax.fill(x_lon,y_lat, color_ton[id])
        if print_id != False:
            x0 = np.mean(x_lon)
            y0 = np.mean(y_lat)
            plt.text(x0, y0, id, fontsize=10)
    if (x_lim != None) & (y_lim != None):     
        plt.xlim(x_lim)
        plt.ylim(y_lim)
    # cax, _ = cbar.make_axes(ax) 
    cax = fig.add_axes([0.45, 0.2, 0.25, 0.03])
    normal = plt.Normalize(data.min(), data.max())
    cb2 = cbar.ColorbarBase(cax, cmap=ListedColormap(color_sq),norm=normal, orientation='horizontal') 
    return ax


def plot_map_fill_multiples_ids(title, city, df,
                                x_lim=None,
                                y_lim=None,
                                figsize=(11, 9),
                                color='r'):
    plt.figure(figsize=figsize)
    fig, ax = plt.subplots(figsize=figsize)
    fig.suptitle(title, fontsize=16)
    for coord in df.coords:
        x = [i[0] for i in coord]
        y = [i[1] for i in coord]
        ax.plot(x, y, 'k')

    for id in city:
        points = df.coords[id]
        x_lon = np.zeros((len(points), 1))
        y_lat = np.zeros((len(points), 1))
        for ip in range(len(points)):
            x_lon[ip] = points[ip][0]
            y_lat[ip] = points[ip][1]
        ax.fill(x_lon, y_lat, color)

        x0 = np.mean(x_lon)
        y0 = np.mean(y_lat)
        plt.text(x0, y0, id, fontsize=10)

    if (x_lim != None) & (y_lim != None):
        plt.xlim(x_lim)
        plt.ylim(y_lim)

    return ax

def plot_shape(sf, id, s=None):
    plt.figure()
    #plotting the graphical axes where map ploting will be done
    ax = plt.axes()
    ax.set_aspect('equal')
#storing the id number to be worked upon
    shape_ex = sf.shape(id)
#NP.ZERO initializes an array of rows and column with 0 in place of each elements
    #an array will be generated where number of rows will be(len(shape_ex,point))and number of columns will be 1 and stored into the variable
    x_lon = np.zeros((len(shape_ex.points),1))
#an array will be generated where number of rows will be(len(shape_ex,point))and number of columns will be 1 and stored into the variable
    y_lat = np.zeros((len(shape_ex.points),1))
    for ip in range(len(shape_ex.points)):
        x_lon[ip] = shape_ex.points[ip][0]
        y_lat[ip] = shape_ex.points[ip][1]
#plotting using the derived coordinated stored in array created by numpy
    plt.plot(x_lon,y_lat)
    x0 = np.mean(x_lon)
    y0 = np.mean(y_lat)
    plt.text(x0, y0, s, fontsize=10)
# use bbox (bounding box) to set plot limits
    plt.xlim(shape_ex.bbox[0],shape_ex.bbox[2])
    return x0, y0

def f_for_plot(f, stpaul_id, tract_sorted, census_tract): # This changes the len of MRT_j from 71 to 82
    # 328 + 329 = 428
    # 348 + 362 = 430
    # 354 + 356 = 429
    # 374[0] = 9800
    for i in range(len(stpaul_id)):
        if census_tract[i] != tract_sorted[i]:
            if tract_sorted[i] in [328, 329]:
                census_tract = np.insert(census_tract, i, 0)
                f = np.insert(f, i, f[-4])
            elif tract_sorted[i] in [348, 362]:
                census_tract = np.insert(census_tract, i, 0)
                f = np.insert(f, i, f[-2])
            elif tract_sorted[i] in [354, 356]:
                census_tract = np.insert(census_tract, i, 0)
                f = np.insert(f, i, f[-3])
            else:
                census_tract = np.insert(census_tract, i, 0)
                f = np.insert(f, i, f[i-1])
    f[-10] = f[-1]
    f = f[:-4]
    return f

def parameter(File_Name='0_2', Lambda=None, Mu=None, e=None, f=None, Pre_L_1=None, Pre_L_2=None, sort=True):
    '''
		:File_Name: File that stored the parameter for the runs
		:Lambda, Mu, e, f, Pre_L_1, Pre_L_2: the date of this collection
		Output: Data
	'''
    file = open(File_Name + '.pkl', 'rb')
    Data = pickle.load(file)
    if Lambda != None:
        Data['Lambda_1'], Data['Lambda_2'] = Lambda
    if Mu != None:
        Data['Mu_1'], Data['Mu_2'] = Mu
    if e != None:
        Data['e'] = np.array(e)
    if f != None:
        Data['f'] = np.array(f)
    if np.any(Pre_L_1 != None):
        Data['Pre_L_1'] = Pre_L_1
    if np.any(Pre_L_2 != None):
        Data['Pre_L_2'] = Pre_L_2
    N, N_1, N_2 = Data['N'], Data['N_1'], Data['N_2']
    if sort:  # If want to perform argsort
        Data['Pre_L_1'] = Data['Pre_L_1'].argsort().argsort()
        Data['Pre_L_2'] = Data['Pre_L_2'].argsort().argsort()
    file.close()
    return Data

def Data_to_Param(Data):
    N_1, N_2 = Data['N_1'], Data['N_2']
    K = Data['K']
    Lambda_1, Lambda_2 = Data['Lambda_1'], Data['Lambda_2']
    Mu_1, Mu_2 = Data['Mu_1'], Data['Mu_2']
    N = Data['N']
    pre_list_1, pre_list_2 = Data['Pre_L_1'], Data['Pre_L_2']
    frac_j_1, frac_j_2 = Data['e'], Data['f']
    return N_1, N_2, K, Lambda_1, Lambda_2, Mu_1, Mu_2, N, pre_list_1, pre_list_2, frac_j_1, frac_j_2

def showPlotFrac(f, x, y, bg, label, bins, text=[], bgcolor=1, color='orange', marker='o', markersize=16):
    # Initializing Visualization Set
    sns.set(style='whitegrid', palette='pastel', color_codes=True)
    sns.mpl.rc('figure', figsize=(10,6))
    #opening the vector map
    shp_path = 'map/Minneapolis_Census_Tracts/Minneapolis_Census_Tracts_wjd.shp'
    # reading the shape file by using reader function of the shape lib
    sf = shp.Reader(shp_path)

    # Converting Shapefile Data Into Pandas Dataframes
    df = read_shapefile(sf)

    tract = df['TRACT'].apply(lambda x:int(x[1:4]))
    stpaul_id = np.where(np.logical_and(tract > 300,tract < 400))[0]
    sorted_stpaul_id = np.argsort(tract[stpaul_id])
    stpaul_df = df.iloc[stpaul_id[sorted_stpaul_id],].reset_index()

    station_y = np.array([44.9517095,44.9772397,44.9760414,44.9436051,44.9655918,44.9517095,44.9772397,44.9760414,44.9270148,44.9567936,44.9332616,44.9446573,44.9781178,44.9560526,44.9037886,44.9705087,44.9493677])
    station_x = np.array([-93.0959001,-93.0320542,-93.1814471,-93.1366808,-93.0566164,-93.0959001,-93.0320542,-93.1814471,-93.1274161,-93.0789894,-93.0827847,-93.1674044,-93.0732163,-93.1291568,-93.1766715,-93.1093789,-93.0256709])
    candidate_x = np.hstack((station_x[:5], station_x[8:]))
    candidate_y = np.hstack((station_y[:5], station_y[8:]))
    text = ['S%s' % (x+1) for x in range(len(candidate_x))]
    title = 'ems'
    tract_sorted = np.sort(tract[stpaul_id])
    census_tract = np.array([301,302,303,304,305,306,307,308,309,310,311,312,313,314,315,316,317,318,319,320,321,322,323,324,325,326,327,330,331,332,333,334,335,336,337,338,339,340,342,344,345,346,347,349,350,351,352,353,355,357,358,359,360,361,363,364,365,366,367,368,369,370,371,372,374,375,376,428,429,430,9800])
    ids = np.arange(len(stpaul_df))
    f_plot = f_for_plot(f, stpaul_id, tract_sorted, census_tract)
    ax = plot_cities_data(stpaul_df, title, ids, f_plot, bgcolor, 15, False, bins)
    ax.plot(station_x, station_y, 'o', marker=marker, color='grey', markersize=markersize, alpha=0.5, label='Candidate Location')
    ax.plot(x, y, 'o', marker=marker, color=color, markersize=markersize, label='Chosen Location')
    for i in range(len(candidate_x)):
        ax.text(candidate_x[i], candidate_y[i], text[i], ha='center', va='center', fontsize=10, weight='bold')

    # if len(text) == len(x):
    #     for i in range(len(x)):
    #         ax.text(x[i], y[i], text[i], ha='center', va='center', fontsize=10)


    # Remove the spines (borders)
    # ax.spines['top'].set_visible(False)
    # ax.spines['right'].set_visible(False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    # Optionally remove ticks
    ax.set_xticks([])



    ax.set_yticks([])

    ax.grid(b=None)
    ax.legend(bbox_to_anchor=(0.627, 0.27), labelspacing=0.8)
    plt.savefig('map_Location_%s_%s.png' %(bg, label), transparent=True, bbox_inches='tight', dpi=400)
    plt.show()

def showPlotIndica(f, x, y, label='Current'):
    # Initializing Visualization Set
    sns.set(style='whitegrid', palette='pastel', color_codes=True)
    sns.mpl.rc('figure', figsize=(10,6))
    #opening the vector map
    shp_path = 'map/Minneapolis_Census_Tracts/Minneapolis_Census_Tracts_wjd.shp'
    # reading the shape file by using reader function of the shape lib
    sf = shp.Reader(shp_path)

    # Converting Shapefile Data Into Pandas Dataframes
    df = read_shapefile(sf)

    tract = df['TRACT'].apply(lambda x:int(x[1:4]))
    stpaul_id = np.where(np.logical_and(tract > 300,tract < 400))[0]
    sorted_stpaul_id = np.argsort(tract[stpaul_id])
    stpaul_df = df.iloc[stpaul_id[sorted_stpaul_id],].reset_index()
    station_y = np.array(
        [44.9517095, 44.9772397, 44.9760414, 44.9436051, 44.9655918, 44.9517095, 44.9772397, 44.9760414, 44.9270148,
         44.9567936, 44.9332616, 44.9446573, 44.9781178, 44.9560526, 44.9037886, 44.9705087, 44.9493677])
    station_x = np.array(
        [-93.0959001, -93.0320542, -93.1814471, -93.1366808, -93.0566164, -93.0959001, -93.0320542, -93.1814471,
         -93.1274161, -93.0789894, -93.0827847, -93.1674044, -93.0732163, -93.1291568, -93.1766715, -93.1093789,
         -93.0256709])

    distance = pd.read_csv("distance.csv")
    distance = distance.drop(distance.columns[0], axis=1)
    distance = distance/1600
    distance[distance > 0.7712735431536174] = (distance[distance > 0.7712735431536174]*111.51113889304331+86.005591195132666)/60
    distance[distance < 0.7712735431536174] = 195.86302790816589*np.sqrt(distance[distance < 0.7712735431536174])/60
    tract_sorted = np.sort(tract[stpaul_id])
    census_tract = np.array(
        [301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322,
         323, 324, 325, 326, 327, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 342, 344, 345, 346, 347, 349,
         350, 351, 352, 353, 355, 357, 358, 359, 360, 361, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 374, 375,
         376, 428, 429, 430, 9800])
    title = 'ems'
    f_plot = f_for_plot(f, stpaul_id, tract_sorted, census_tract)
    cities = np.where(np.array(f_plot) == 1)[0]
    ax = plot_map_fill_multiples_ids('', cities, stpaul_df, color = 'g')
    # ax = plot_cities_data(stpaul_df, title, ids, f_plot, 1, 15, False)
    ax.plot(station_x[3:],station_y[3:], color='orange', marker='o', markersize=12, ls='', label = 'Candidate')
    ax.plot(x, y, 'bs', markersize=12, label = label)
    ax.legend()
    ax.grid(b=None)
    plt.savefig('map_Location_indi_%s.png' %label, transparent=True,bbox_inches='tight', dpi=400)
    plt.show()

if __name__ == '__main__':

    data = pd.read_pickle('C:\\Users\\Yixing\\Dropbox\\Research\\SpatialQueue\\(2nd)POM\\OperDecision\\data_9units.pkl')
    Data = parameter(File_Name='Saint_Paul')
    N_1, N_2, K, Lambda_1, Lambda_2, Mu_1, Mu_2, N, pre_list_1, pre_list_2, frac_j_1, frac_j_2 = Data_to_Param(Data)
    f = frac_j_1
    station_y = np.array([44.9517095,44.9772397,44.9760414,44.9436051,44.9655918,44.9517095,44.9772397,44.9760414,44.9270148,44.9567936,44.9332616,44.9446573,44.9781178,44.9560526,44.9037886,44.9705087,44.9493677])
    station_x = np.array([-93.0959001,-93.0320542,-93.1814471,-93.1366808,-93.0566164,-93.0959001,-93.0320542,-93.1814471,-93.1274161,-93.0789894,-93.0827847,-93.1674044,-93.0732163,-93.1291568,-93.1766715,-93.1093789,-93.0256709])
    # candidate_x = np.hstack((station_x[:5], station_x[8:]))
    # candidate_y = np.hstack((station_y[:5], station_y[8:]))
    # stations = ['S%s' %x for x in range(14)]
    # showPlotFrac(f, candidate_x, candidate_y, text=stations, bg='fj', label='Candidate', bgcolor=1, marker='s', color='b')

    currentLoc = data.loc[0]
    currentNo = currentLoc['location'].copy()
    currentNo = [x - 3 if x > 5 else x for x in currentNo]
    currentNo = ['S%s' %(x+1) for x in currentNo]
    currentLoc_x = [station_x[i] for i in currentLoc['location']]
    currentLoc_y = [station_y[i] for i in currentLoc['location']]
    currentDist_F6 = currentLoc['frac_6'][:-1].copy()
    currentDist_F9 = currentLoc['frac_9'][:-1].copy()
    currentMRT = currentLoc['MRT_j_9'].copy()
    currentMRT_f = currentLoc['MRT_j_9'].copy() * f

    bestLoc_6 = data.loc[data['frac_total_4.5'].idxmax()]
    bestNo_6 = bestLoc_6['location'].copy()
    bestNo_6 = [x - 3 if x > 5 else x for x in bestNo_6]
    bestNo_6 = ['S%s' % (x + 1) for x in bestNo_6]
    bestLoc_6_x = [station_x[i] for i in bestLoc_6['location']]
    bestLoc_6_y = [station_y[i] for i in bestLoc_6['location']]
    bestDist_6_F6 = bestLoc_6['frac_6'][:-1].copy()
    bestDist_6_F9 = bestLoc_6['frac_9'][:-1].copy()
    bestMRT_6 = bestLoc_6['MRT_j_6'].copy()
    bestMRT_6_f = bestLoc_6['MRT_j_6'].copy() * f

    bestLoc_9 = data.loc[data['frac_total_7.5'].idxmax()]
    bestNo_9 = bestLoc_9['location'].copy()
    bestNo_9 = [x - 3 if x > 5 else x for x in bestNo_9]
    bestNo_9 = ['S%s' % (x+1) for x in bestNo_9]
    bestLoc_9_x = [station_x[i] for i in bestLoc_9['location']]
    bestLoc_9_y = [station_y[i] for i in bestLoc_9['location']]
    bestDist_9_F6 = bestLoc_9['frac_6'][:-1].copy()
    bestDist_9_F9 = bestLoc_9['frac_9'][:-1].copy()
    bestMRT_9 = bestLoc_9['MRT_j_9'].copy()
    bestMRT_9_f = bestLoc_9['MRT_j_9'].copy() * f

    bestLoc_MRT = data.loc[data['MRT'].idxmin()]
    bestNo_MRT = bestLoc_MRT['location'].copy()
    bestNo_MRT = [x - 3 if x > 5 else x for x in bestNo_MRT]
    bestNo_MRT = ['S%s' % (x+1) for x in bestNo_MRT]
    bestLoc_MRT_x = [station_x[i] for i in bestLoc_MRT['location']]
    bestLoc_MRT_y = [station_y[i] for i in bestLoc_MRT['location']]
    bestDist_MRT_F6 = bestLoc_MRT['frac_6'][:-1].copy()
    bestDist_MRT_F9 = bestLoc_MRT['frac_9'][:-1].copy()
    bestMRT_MRT = bestLoc_MRT['MRT_j_9'].copy()
    bestMRT_MRT_f = bestLoc_MRT['MRT_j_9'].copy() * f
    #
    # # # Calculate the bins
    # MRT_bins = np.hstack((currentMRT, bestMRT_MRT))
    # new_MRT_bins, bins_MRT = pd.qcut(MRT_bins, 15, retbins=True, labels=False, duplicates='drop')
    # np.save('bins_MRT.npy', bins_MRT)
    # showPlotFrac(currentMRT, currentLoc_x, currentLoc_y, bins=bins_MRT, text=currentNo, bg='MRT', label='Current', bgcolor=4)
    # showPlotFrac(bestMRT_MRT, bestLoc_MRT_x, bestLoc_MRT_y, bins=bins_MRT, text=bestNo_MRT, bg='MRT', label='Best', bgcolor=4)

    # MRTj_bins = np.hstack((currentMRT_f, bestMRT_MRT_f))
    # new_MRTj_bins, bins_MRT_f = pd.qcut(MRTj_bins, 15, retbins=True, labels=False, duplicates='drop')
    # np.save('bins_MRT_f.npy', bins_MRT_f)
    # showPlotFrac(currentMRT_f, currentLoc_x, currentLoc_y, bins=bins_MRT_f, text=currentNo, bg='Weighted MRT',
    #              label='Current', bgcolor=4)
    # showPlotFrac(bestMRT_MRT_f, bestLoc_MRT_x, bestLoc_MRT_y, bins=bins_MRT_f, text=bestNo_MRT, bg='Weighted MRT', label='Best', bgcolor=4)

    # # # Calculate the bins
    # f6_bins = np.hstack((currentDist_F6, bestDist_6_F6))
    # new_f6_bins, bins_F6 = pd.qcut(f6_bins, 15, retbins=True, labels=False, duplicates='drop')
    # np.save('bins_F6.npy', bins_F6)
    # showPlotFrac(currentDist_F6, currentLoc_x, currentLoc_y, bins=bins_F6, text=currentNo, bg='F(6)', label='Current', bgcolor=2)
    # showPlotFrac(bestDist_6_F6, bestLoc_6_x, bestLoc_6_y, bins=bins_F6, text=bestNo_6, bg='F(6)', label='Best', bgcolor=2)
    #
    f9_bins = np.hstack((currentDist_F9, bestDist_9_F9))
    new_f9_bins, bins_F9 = pd.qcut(f9_bins, 15, retbins=True, labels=False, duplicates='drop')
    np.save('bins_F9.npy', bins_F9)
    showPlotFrac(currentDist_F9, currentLoc_x, currentLoc_y, bins=bins_F9, text=currentNo, bg='F(9)', label='Current', bgcolor=8)
    showPlotFrac(bestDist_9_F9, bestLoc_9_x, bestLoc_9_y, bins=bins_F9, text=bestNo_9, bg='F(9)', label='Best', bgcolor=8)

    # fig, ax = plt.subplots()
    # x = np.arange(len(currentMRT)) + 1
    # y = list(currentMRT - bestMRT_6)
    # y.sort(reverse=True)
    # ax.bar(x, y, color='black')
    # ax.set_xlim(0, 72)
    # ax.set_xlabel("Index of Nodes")
    # ax.set_ylabel('Difference of $MRT_j$')
    # yposition = np.arange(-6.0, 5.0, 1.0)
    # xposition = np.arange(1, 71, 5.0)
    # for yc in yposition:
    #     plt.axhline(y=yc, color='lightgrey', linestyle='-.', linewidth=0.7)
    # for xc in xposition:
    #     plt.axvline(x=xc, color='lightgrey', linestyle='-.', linewidth=0.7)
    # # plt.grid(b=False)
    # plt.xticks(np.arange(1, 72, 5), fontsize=10)
    # plt.yticks(np.arange(-6, 5.1, 1), fontsize=10)
    # # plt.legend(loc=1)
    # plt.tight_layout()
    # plt.savefig('bar_MRT_F6_diff.png', transparent=False, bbox_inches='tight', dpi=400)
    # plt.show()
    #
    # fig, ax = plt.subplots()
    # x = np.arange(len(currentMRT_f)) + 1
    # y = list(currentMRT_f - bestMRT_6_f)
    # y.sort(reverse=True)
    # ax.bar(x, y, color='black')
    # ax.set_xlim(0, 72)
    # ax.set_xlabel("Index of Atoms")
    # ax.set_ylabel('Difference of $MRT_j$')
    # yposition = np.arange(-0.2, 0.2, 0.05)
    # xposition = np.arange(1, 71, 5.0)
    # for yc in yposition:
    #     plt.axhline(y=yc, color='lightgrey', linestyle='-.', linewidth=0.7)
    # for xc in xposition:
    #     plt.axvline(x=xc, color='lightgrey', linestyle='-.', linewidth=0.7)
    # # plt.grid(b=False)
    # plt.xticks(np.arange(1, 72, 5), fontsize=10)
    # plt.yticks(np.arange(-0.2, 0.2, 0.05), fontsize=10)
    # # plt.legend(loc=1)
    # plt.tight_layout()
    # plt.savefig('bar_MRT_f_F6_diff.png', transparent=False, bbox_inches='tight', dpi=400)
    # plt.show()
    #
    # fig, ax = plt.subplots()
    # x = np.arange(len(currentDist_F6))
    # y = currentDist_F6 - bestDist_6_F6
    # y.sort()
    # n, bins, patches = ax.hist(y, color='black', bins=np.linspace(-0.9, 0.9, 19))
    # for i in range(len(patches)):
    #     height = patches[i].get_height()
    #     if height > 0:
    #         # Calculate the position for text (center of the bin)
    #         bin_center = 0.5 * (bins[i] + bins[i + 1])
    #         # Place text at the bottom of the bar
    #         ax.text(bin_center, height, f'{int(height)}', ha='center', va='bottom', fontsize=10, color='black')
    # ax.set_xlim(-0.9,0.9)
    # ax.set_ylabel("Number of Atoms")
    # ax.set_xlabel('Difference of $F_j(6)$')
    # xposition = np.arange(-0.9, 0.9, 0.1)
    # for xc in xposition:
    #     plt.axvline(x=xc, color='lightgrey', linestyle='-.', linewidth=0.7)
    # plt.grid(b=False)
    # plt.xticks(np.arange(-0.9, 0.9, 0.1), rotation=60, fontsize=10)
    # # plt.legend(loc=1)
    # plt.tight_layout()
    # plt.savefig('hist_F6_F6_diff.png', transparent=False, bbox_inches='tight', dpi=400)
    # plt.show()

    # fig, ax = plt.subplots()
    # x = np.arange(len(currentDist_F9))
    # y = currentDist_F9 - bestDist_6_F9
    # y.sort()
    # n, bins, patches = ax.hist(y, color='black', bins=np.linspace(-0.9, 0.9, 19))
    # for i in range(len(patches)):
    #     height = patches[i].get_height()
    #     if height > 0:
    #         # Calculate the position for text (center of the bin)
    #         bin_center = 0.5 * (bins[i] + bins[i + 1])
    #         # Place text at the bottom of the bar
    #         ax.text(bin_center, height, f'{int(height)}', ha='center', va='bottom', fontsize=10, color='black')
    # ax.set_xlim(-0.9,0.9)
    # ax.set_ylabel("Number of Atoms")
    # ax.set_xlabel('Difference of $F_j(9)$')
    # xposition = np.arange(-0.9, 0.9, 0.1)
    # for xc in xposition:
    #     plt.axvline(x=xc, color='lightgrey', linestyle='-.', linewidth=0.7)
    # plt.grid(b=False)
    # plt.xticks(np.arange(-0.9, 0.9, 0.1), rotation=60, fontsize=10)
    # # plt.legend(loc=1)
    # plt.tight_layout()
    # plt.savefig('hist_F9_F6_diff.png', transparent=False, bbox_inches='tight', dpi=400)
    # plt.show()

    # fig, ax = plt.subplots()
    # x = np.arange(len(currentMRT)) + 1
    # y = list(currentMRT - bestMRT_9)
    # y.sort(reverse=True)
    # ax.bar(x, y, color='black')
    # ax.set_xlim(0, 72)
    # ax.set_xlabel("Index of Atoms")
    # ax.set_ylabel('Difference of $MRT_j$')
    # yposition = np.arange(-3, 3, 0.5)
    # xposition = np.arange(1, 71, 5.0)
    # for yc in yposition:
    #     plt.axhline(y=yc, color='lightgrey', linestyle='-.', linewidth=0.7)
    # for xc in xposition:
    #     plt.axvline(x=xc, color='lightgrey', linestyle='-.', linewidth=0.7)
    # # plt.grid(b=False)
    # plt.xticks(np.arange(1, 72, 5), fontsize=10)
    # plt.yticks(np.arange(-3, 3, 0.5), fontsize=10)
    # # plt.legend(loc=1)
    # plt.tight_layout()
    # plt.savefig('bar_MRT_F9_diff.png', transparent=False, bbox_inches='tight', dpi=400)
    # plt.show()
    #
    # fig, ax = plt.subplots()
    # x = np.arange(len(currentMRT_f)) + 1
    # y = list(currentMRT_f - bestMRT_9_f)
    # y.sort(reverse=True)
    # ax.bar(x, y, color='black')
    # ax.set_xlim(0, 72)
    # ax.set_xlabel("Index of Atoms")
    # ax.set_ylabel('Difference of $MRT_j$')
    # yposition = np.arange(-0.06, 0.06, 0.02)
    # xposition = np.arange(1, 71, 5.0)
    # for yc in yposition:
    #     plt.axhline(y=yc, color='lightgrey', linestyle='-.', linewidth=0.7)
    # for xc in xposition:
    #     plt.axvline(x=xc, color='lightgrey', linestyle='-.', linewidth=0.7)
    # # plt.grid(b=False)
    # plt.xticks(np.arange(1, 72, 5), fontsize=10)
    # plt.yticks(np.arange(-0.06, 0.06, 0.02), fontsize=10)
    # # plt.legend(loc=1)
    # plt.tight_layout()
    # plt.savefig('bar_MRT_f_F9_diff.png', transparent=False, bbox_inches='tight', dpi=400)
    # plt.show()

    # fig, ax = plt.subplots()
    # x = np.arange(len(currentDist_F6))
    # y = currentDist_F6 - bestDist_9_F6
    # y.sort()
    # n, bins, patches = ax.hist(y, color='black', bins=np.linspace(-0.9, 0.9, 19))
    # for i in range(len(patches)):
    #     height = patches[i].get_height()
    #     if height > 0:
    #         # Calculate the position for text (center of the bin)
    #         bin_center = 0.5 * (bins[i] + bins[i + 1])
    #         # Place text at the bottom of the bar
    #         ax.text(bin_center, height, f'{int(height)}', ha='center', va='bottom', fontsize=10, color='black')
    # ax.set_xlim(-0.9,0.9)
    # ax.set_ylabel("Number of Atoms")
    # ax.set_xlabel('Difference of $F_j(6)$')
    # xposition = np.arange(-0.9, 0.9, 0.1)
    # for xc in xposition:
    #     plt.axvline(x=xc, color='lightgrey', linestyle='-.', linewidth=0.7)
    # plt.grid(b=False)
    # plt.xticks(np.arange(-0.9, 0.9, 0.1), rotation=60, fontsize=10)
    # # plt.legend(loc=1)
    # plt.tight_layout()
    # plt.savefig('hist_F6_F9_diff.png', transparent=False, bbox_inches='tight', dpi=400)
    # plt.show()

    # fig, ax = plt.subplots()
    # x = np.arange(len(currentDist_F9))
    # y = currentDist_F9 - bestDist_9_F9
    # y.sort()
    # n, bins, patches = ax.hist(y, color='black', bins=np.linspace(-0.9, 0.9, 19))
    # for i in range(len(patches)):
    #     height = patches[i].get_height()
    #     if height > 0:
    #         # Calculate the position for text (center of the bin)
    #         bin_center = 0.5 * (bins[i] + bins[i + 1])
    #         # Place text at the bottom of the bar
    #         ax.text(bin_center, height, f'{int(height)}', ha='center', va='bottom', fontsize=10, color='black')
    # ax.set_xlim(-0.9,0.9)
    # ax.set_ylabel("Number of Atoms")
    # ax.set_xlabel('Difference of $F_j(9)$')
    # xposition = np.arange(-0.9, 0.9, 0.1)
    # for xc in xposition:
    #     plt.axvline(x=xc, color='lightgrey', linestyle='-.', linewidth=0.7)
    # plt.grid(b=False)
    # plt.xticks(np.arange(-0.9, 0.9, 0.1), rotation=60, fontsize=10)
    # # plt.legend(loc=1)
    # plt.tight_layout()
    # plt.savefig('hist_F9_F9_diff.png', transparent=False, bbox_inches='tight', dpi=400)
    # plt.show()

    # fig, ax = plt.subplots()
    # x = np.arange(len(currentMRT)) + 1
    # y = list(currentMRT - bestMRT_MRT)
    # y.sort(reverse=True)
    # ax.bar(x, y, color='black')
    # ax.set_xlim(0, 72)
    # ax.set_xlabel("Index of Nodes")
    # ax.set_ylabel('Difference of $MRT_j$')
    # yposition = np.arange(-3.0, 3.0, 1.0)
    # xposition = np.arange(1, 71, 5.0)
    # for yc in yposition:
    #     plt.axhline(y=yc, color='lightgrey', linestyle='-.', linewidth=0.7)
    # for xc in xposition:
    #     plt.axvline(x=xc, color='lightgrey', linestyle='-.', linewidth=0.7)
    # # plt.grid(b=False)
    # plt.xticks(np.arange(1, 72, 5), fontsize=10)
    # plt.yticks(np.arange(-3.0, 3.1, 1), fontsize=10)
    # # plt.legend(loc=1)
    # plt.tight_layout()
    # plt.savefig('bar_MRT_MRT_diff.png', transparent=False, bbox_inches='tight', dpi=400)
    # plt.show()
    # # #
    # fig, ax = plt.subplots()
    # x = np.arange(len(currentMRT_f)) + 1
    # y = list(currentMRT_f - bestMRT_MRT_f)
    # y.sort(reverse=True)
    # ax.bar(x, y, color='black')
    # ax.set_xlim(0, 72)
    # ax.set_xlabel("Index of Atoms")
    # ax.set_ylabel('Difference of Weighted $MRT_j$')
    # yposition = np.arange(-0.06, 0.06, 0.01)
    # xposition = np.arange(1, 71, 5.0)
    # for yc in yposition:
    #     plt.axhline(y=yc, color='lightgrey', linestyle='-.', linewidth=0.7)
    # for xc in xposition:
    #     plt.axvline(x=xc, color='lightgrey', linestyle='-.', linewidth=0.7)
    # # plt.grid(b=False)
    # plt.xticks(np.arange(1, 72, 5), fontsize=10)
    # plt.yticks(np.arange(-0.06, 0.06, 0.01), fontsize=10)
    # # plt.legend(loc=1)
    # plt.tight_layout()
    # plt.savefig('bar_MRT_f_MRT_diff.png', transparent=False, bbox_inches='tight', dpi=400)
    # plt.show()
    #
    # fig, ax = plt.subplots()
    # x = np.arange(len(currentDist_F6))
    # y = currentDist_F6 - bestDist_MRT_F6
    # y.sort()
    # n, bins, patches = ax.hist(y, color='black', bins=np.linspace(-0.9, 0.9, 19))
    # for i in range(len(patches)):
    #     height = patches[i].get_height()
    #     if height > 0:
    #         # Calculate the position for text (center of the bin)
    #         bin_center = 0.5 * (bins[i] + bins[i + 1])
    #         # Place text at the bottom of the bar
    #         ax.text(bin_center, height, f'{int(height)}', ha='center', va='bottom', fontsize=10, color='black')
    # ax.set_xlim(-0.9,0.9)
    # ax.set_ylabel("Number of Atoms")
    # ax.set_xlabel('Difference of $F_j(6)$')
    # xposition = np.arange(-0.9, 0.9, 0.1)
    # for xc in xposition:
    #     plt.axvline(x=xc, color='lightgrey', linestyle='-.', linewidth=0.7)
    # plt.grid(b=False)
    # plt.xticks(np.arange(-0.9, 0.9, 0.1), rotation=60, fontsize=10)
    # # plt.legend(loc=1)
    # plt.tight_layout()
    # plt.savefig('hist_F6_MRT_diff.png', transparent=False, bbox_inches='tight', dpi=400)
    # plt.show()
    #
    # fig, ax = plt.subplots()
    # x = np.arange(len(currentDist_F9))
    # y = currentDist_F9 - bestDist_MRT_F9
    # y.sort()
    # n, bins, patches = ax.hist(y, color='black', bins=np.linspace(-0.9, 0.9, 19))
    # for i in range(len(patches)):
    #     height = patches[i].get_height()
    #     if height > 0:
    #         # Calculate the position for text (center of the bin)
    #         bin_center = 0.5 * (bins[i] + bins[i + 1])
    #         # Place text at the bottom of the bar
    #         ax.text(bin_center, height, f'{int(height)}', ha='center', va='bottom', fontsize=10, color='black')
    # ax.set_xlim(-0.9,0.9)
    # ax.set_ylabel("Number of Atoms")
    # ax.set_xlabel('Difference of $F_j(9)$')
    # xposition = np.arange(-0.9, 0.9, 0.1)
    # for xc in xposition:
    #     plt.axvline(x=xc, color='lightgrey', linestyle='-.', linewidth=0.7)
    # plt.grid(b=False)
    # plt.xticks(np.arange(-0.9, 0.9, 0.1), rotation=60, fontsize=10)
    # # plt.legend(loc=1)
    # plt.tight_layout()
    # plt.savefig('hist_F9_MRT_diff.png', transparent=False, bbox_inches='tight', dpi=400)
    # plt.show()