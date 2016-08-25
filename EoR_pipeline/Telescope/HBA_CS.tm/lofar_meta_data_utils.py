
def get_ant_positions(file):
    import numpy as np;
    import re;
    lba=None;lba_rot=None;hba=None;hba0=None;hba1=None;hba_rot=None;hba0_rot=None;hba1_rot=None;lba_centre=None;hba_centre=None;hba0_centre=None;hba1_centre=None;
    fl=open(file);
    fl.seek(0);
    lba_rot_re = re.compile("ROTATION_MATRIX LBA");
    lines=fl.readlines();
    for iline,line in enumerate(lines):
        if lba_rot_re.match(line):
            matrix_size = np.array(re.findall("\d+",lines[iline+1]),np.int32);   
            lba_rot =np.zeros(matrix_size);
            for i in range(matrix_size[0]):
                lba_rot[i]=np.array(re.findall("-?\d+\.\d+",lines[iline+2+i]),np.float).reshape(matrix_size[1:]);
            break;    
    fl.seek(0);
    lba_re = re.compile("LBA\n");
    for iline,line in enumerate(lines):
        if lba_re.match(line):
            lba_centre=np.array(re.findall("\d+\.\d+",lines[iline+1]),np.float);
            matrix_size = np.array(re.findall("\d+",lines[iline+2]),np.int32);   
            lba = np.zeros(matrix_size);
            for i in range(matrix_size[0]):
                a=re.findall("-?\d+\.\d+",lines[iline+3+i]);
                lba[i]=np.array(a).reshape(matrix_size[1:]);
            lba+=lba_centre;
            break;
    hba_re = re.compile("HBA\n");
    fl.seek(0);
    for iline,line in enumerate(lines):
        if hba_re.match(line):
            hba_centre=np.array(re.findall("\d+\.\d+",lines[iline+1]),np.float);
            matrix_size = np.array(re.findall("\d+",lines[iline+2]),np.int32);   
            matrix_size[0]/=2;
            hba0 = np.zeros(matrix_size);
            hba1 = np.zeros(matrix_size);
            for i in range(matrix_size[0]):
                a=re.findall("-?\d+\.\d+",lines[iline+3+i]);
                hba0[i]=np.array(a).reshape(matrix_size[1:]);
            for i in range(matrix_size[0],2*matrix_size[0]):
                a=re.findall("-?\d+\.\d+",lines[iline+3+i]);
                hba1[i-matrix_size[0]]=np.array(a).reshape(matrix_size[1:]);
            hba0+=hba_centre;
            hba1+=hba_centre;
            hba = np.vstack((hba0,hba1));
            break;
    hba_rot_re = re.compile("ROTATION_MATRIX HBA0");
    fl.seek(0);
    for iline,line in enumerate(lines):
        if hba_rot_re.match(line):
            matrix_size = np.array(re.findall("\d+",lines[iline+1]),np.int32);   
            hba0_rot =np.zeros(matrix_size);
            for i in range(matrix_size[0]):
                hba0_rot[i]=np.array(re.findall("-?\d+\.\d+",lines[iline+2+i]),np.float).reshape(matrix_size[1:]);
            break;   
    hba_rot_re = re.compile("ROTATION_MATRIX HBA1");
    fl.seek(0);
    for iline,line in enumerate(lines):
        if hba_rot_re.match(line):
            matrix_size = np.array(re.findall("\d+",lines[iline+1]),np.int32);   
            hba1_rot =np.zeros(matrix_size);
            for i in range(matrix_size[0]):
                hba1_rot[i]=np.array(re.findall("-?\d+\.\d+",lines[iline+2+i]),np.float).reshape(matrix_size[1:]);
            break;
    hba_rot_re = re.compile("ROTATION_MATRIX HBA");
    fl.seek(0);
    for iline,line in enumerate(lines):
        if hba_rot_re.match(line):
            matrix_size = np.array(re.findall("\d+",lines[iline+1]),np.int32);   
            hba_rot =np.zeros(matrix_size);
            for i in range(matrix_size[0]):
                hba_rot[i]=np.array(re.findall("-?\d+\.\d+",lines[iline+2+i]),np.float).reshape(matrix_size[1:]);
            break;   
    hba_re = re.compile("HBA0\n");
    fl.seek(0);
    for iline,line in enumerate(lines):
        if hba_re.match(line):
            hba0_centre=np.array(re.findall("\d+\.\d+",lines[iline+1]),np.float);
            break;
    hba_re = re.compile("HBA1\n");
    fl.seek(0);
    for iline,line in enumerate(lines):
        if hba_re.match(line):
            hba1_centre=np.array(re.findall("\d+\.\d+",lines[iline+1]),np.float);
            break;
    fl.close();
    return (lba,lba_rot,hba,hba0,hba1,hba_rot,hba0_rot,hba1_rot,lba_centre,hba_centre,hba0_centre,hba1_centre);


def get_hba_deltas(file):
    import numpy as np;
    import re;
    fl=open(file);
    hba_deltas_re = re.compile("HBADeltas");
    lines = fl.readlines();
    for iline,line in enumerate(lines):
        if hba_deltas_re.match(line):
            matrix_size = np.array(re.findall("\d+",lines[iline+1]),np.int32);   
            hba_deltas =np.zeros(matrix_size);
            for i in range(matrix_size[0]):
                hba_deltas[i]=np.array(re.findall("-?\d+\.\d+",lines[iline+2+i]),np.float).reshape(matrix_size[1:]);
    fl.close();
    return hba_deltas;


def load_data(path="./StaticMetaData/"):
    import os;
    import re;
    filelist = os.listdir(path);
    if os.access(path+"AntennaFields/",os.F_OK):
        filelist +=  os.listdir(path+"AntennaFields/");
    if os.access(path+"iHBADeltas/",os.F_OK):
        filelist +=  os.listdir(path+"iHBADeltas/");
    ant_conf = re.compile("(.{5})(-AntennaField\.conf)");
    hba_conf = re.compile("(.{5})(-iHBADeltas\.conf)");
    ant_list = [];
    hba_list = [];
    station_dict = [];
    hba_deltas = {};
    station_dict = {};
    for file in filelist:
        am = ant_conf.match(file);
        if am:
            ant_list.append(am.group(1));
            #print "getting data for station",am.group(1);
            if os.access(path+"AntennaFields/"+file,os.F_OK):
                (lba,lba_rot,hba,hba0,hba1,hba_rot,hba0_rot,hba1_rot,lba_c,hba_c,hba0_c,hba1_c) = get_ant_positions(path+"AntennaFields/"+file);
            else:
                (lba,lba_rot,hba,hba0,hba1,hba_rot,hba0_rot,hba1_rot,lba_c,hba_c,hba0_c,hba1_c) = get_ant_positions(path+file);
            station_dict[am.group(1)]=dict(lba=lba,lba_c_rot=lba_rot,hba=hba,hba0=hba0,hba1=hba1,hba_c_rot=hba_rot,hba0_c_rot=hba0_rot,hba1_c_rot=hba1_rot,lba_c=lba_c,hba_c=hba_c,hba0_c=hba0_c,hba1_c=hba1_c);
        else:
            hm = hba_conf.match(file);
            if hm:
                #print "getting hba deltas for station",hm.group(1),path+"iHBADeltas/"+file;
                hba_list.append(hm.group(1));
                if os.access(path+"iHBADeltas/"+file,os.F_OK):
                    hba_deltas[hm.group(1)]=get_hba_deltas(path+"iHBADeltas/"+file);
                else:
                    hba_deltas[hm.group(1)]=get_hba_deltas(path+file);
    ant_list = set(ant_list);
    hba_list= set(hba_list)&ant_list; #stations that have hba information
    for station in hba_list:
        station_dict[station]['hba_deltas']=hba_deltas[station];
    return station_dict;


