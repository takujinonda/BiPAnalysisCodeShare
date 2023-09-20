# -*- coding: utf-8 -*-
import argparse
import pandas as pd
from main_func import main_func

def get_args():
    parser = argparse.ArgumentParser(
            prog='BiP Analysis',
            usage='Input url of csv data',
            description='Input url of csv data',
            epilog='end',
            add_help=True,
            )
    
    parser.add_argument('data_url', help='url of csv data', type = str)
    args = parser.parse_args()

    return(args)


def main():
    args = get_args()

    data_url = args.data_url
    
    AnalyzeSensorName = ['time', 'latitude', 'longitude']

    # datapq = pq.read_table(data_url, filesystem=s3, columns=['latitude', 'longitude', 'time'])
    # datapq = pq.read_table(data_url, columns=AnalyzeSensorName)
    # data = datapq.to_pandas()
    data = pd.read_csv(data_url, usecols=AnalyzeSensorName)

    outdf = main_func(data)
    outdf.to_csv("out.csv")


if __name__ == '__main__':
    main()

