#include "external.h"
#include "stdio.h"

#define SRC_PATH "/home/linzhiqiang/D/buildMapping_CPP/map_R_new_undistort/yuvImage/test"

void main()
{
    unsigned char *image_data = (unsigned char*)malloc((640*480));
    char *img_src = malloc(128);

    int run =1;
    int index = 351;
    FILE *fd;
    char m_bool= 0;
    init_map_src();

    while(run == 1)
    {
        sprintf(img_src,"%s/%05d.yuv",SRC_PATH,index);
        fd = fopen(img_src, "rb");
        int readsize = fread(image_data, 1, 480*640, fd);
        fclose(fd);
        printf("img_src %s  %d\n",img_src,readsize);

        if(index == 1160)
        {
            m_bool = 1;
            run = 0;
        }
        maping_world(480,640,image_data,m_bool);
        index++;
    }
    printf("end !!!!!!!!!!! \n");

    free(image_data);
    free(img_src);

}