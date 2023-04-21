from futures3.thread import ThreadPoolExecutor as Multithread
#from requests import Request, Session
import requests
from requests.exceptions import ConnectionError, Timeout, TooManyRedirects
from pytrends.request import TrendReq
import yfinance as yf
import pandas as pd
import openpyxl
import finviz
import time
import json

MAX_THREADS = 30 
API_KEY_CMC = 'API_KEY_CMC'
API_KEY_CRYPTO_COMPARE = 'API_KEY_CRYPTO_COMPARE'
API_KEY_BINANCE = 'API_KEY_BINANCE'

def flatten_dict(my_dict, existing_dict):
    """Flatten A Given Dictionary
    :param my_dict: dict containing nested dict
    :type my_dict: dict
    :param existing_dict: dict to append flattened data
    :type existing_dict: dict
    :return: dict containing flattened data
    """
    for key, value in my_dict.items():
        if not isinstance(value, dict):
            existing_dict[key] = value
        else:
            flatten_dict(value, existing_dict)
    return existing_dict

def get_value_by_named_range(strPath, strNamedRange):
    """Get Given Named Range Cell From Given Excel File
    :param strPath: Full path of excel file
    :type strPath: string
    :param strNamedRange: Named Range
    :type strNamedRange: string
    :return: Cell Object
    """
    wb = openpyxl.load_workbook(strPath, data_only = True)
    my_range = wb.defined_names[strNamedRange] 
    dest = my_range.destinations
    for title, coord in dest:
        ws = wb[title]
        cells = ws[coord]
        
    return cells, ws

def get_crypto_pricing(tickers_string, strFileName):
    url = 'https://min-api.cryptocompare.com/data/pricemultifull'
    parameters = {
        'fsyms': tickers_string,
        'tsyms':'USD'
    }
    headers = {
    'Accepts': 'application/json',
    'authorization': API_KEY_CRYPTO_COMPARE
    }

    session = requests.Session()
    session.headers.update(headers)

    try:
        response = session.get(url, params = parameters)
        data = json.loads(response.text)
        data_flat = {}
        for key, value in data['DISPLAY'].items():
            data_flat[key] = flatten_dict(value,{})
        crypto_df = pd.DataFrame.from_dict(data_flat,orient = 'index')
        crypto_df.reset_index(inplace = True)
        crypto_df = crypto_df.rename(columns = {'index':'SYMBOL'})
        crypto_df.to_csv(path_or_buf = strFileName,index = False)
    except (ConnectionError, Timeout, TooManyRedirects) as e:
        print(e)

def get_crypto_global_metrics(strFileName):
    url = 'https://pro-api.coinmarketcap.com/v1/global-metrics/quotes/latest'
    parameters = {
    'convert':'USD'
    }
    headers = {
    'Accepts': 'application/json',
    'X-CMC_PRO_API_KEY': API_KEY_CMC,
    }

    session = requests.Session()
    session.headers.update(headers)

    try:
        response = session.get(url, params = parameters)
        data = json.loads(response.text)
        data_flat = {}
        for key, value in data.items():
            data_flat[key] = flatten_dict(value,{})
        crypto_df = pd.DataFrame.from_dict(data_flat,orient = 'index')
        crypto_df.to_csv(path_or_buf = strFileName,index = False)
    except (ConnectionError, Timeout, TooManyRedirects) as e:
        print(e)

def get_goog_trends(keywords_list, strFileName):
    pytrend = TrendReq()
    kw_list = keywords_list
    pytrend.build_payload(kw_list,timeframe = 'today 5-y')
    df_trend = pytrend.interest_over_time()
    df_trend.reset_index(inplace = True)
    df_trend.to_csv(path_or_buf = strFileName,index = False)

def get_tx_fees(strFileName):
    strUrl = 'https://api.blockchair.com/bitcoin/transactions?a=date,avg(fee_usd)&export=csv'
    req = requests.get(strUrl)
    url_content = req.content
    csv_file = open(strFileName, 'wb')
    csv_file.write(url_content)
    csv_file.close()

def get_futures_binance(strFileName):
    Binance_BaseUrl = 'https://fapi.binance.com'
    Binance_Premium_Endpoint = '/fapi/v1/premiumIndex'
    Binance_OI_Endpoint = '/fapi/v1/openInterest'

    parameters = {
        'SYMBOL':'BTCUSDT'
    }
    session = requests.Session()

    try:
        response = session.get(Binance_BaseUrl + Binance_Premium_Endpoint, params = parameters)
        data_0 = json.loads(response.text)
        response = session.get(Binance_BaseUrl + Binance_OI_Endpoint, params = parameters)
        data_1 = json.loads(response.text)
        data_0['openInterest'] = data_1['openInterest']
        df_futures = pd.DataFrame.from_dict(data_0,orient='index').transpose()
        df_futures.to_csv(path_or_buf = strFileName,index = False)
    except (ConnectionError, Timeout, TooManyRedirects) as e:
        print(e)

if __name__ == '__main__':
    #GET Historical and Current Pricing Data
    t0 = time.time()
    path = r'C:\Users\JustinCata.000\PortfolioMonitor\PriceAlertMonitor.xlsm'

    cells, ws = get_value_by_named_range(path, 'CC_TickerList')
    cells_toCSV_FileName, ws = get_value_by_named_range(path, 'CC_FullPath')
    get_crypto_pricing(cells.value, cells_toCSV_FileName.value)

    cells_toCSV_FileName, ws = get_value_by_named_range(path, 'CoinMC_GM_FullPath')
    get_crypto_global_metrics(cells_toCSV_FileName.value)

    cells, ws = get_value_by_named_range(path, 'Goog_WordList')
    cells_toCSV_FileName, ws = get_value_by_named_range(path, 'Goog_FullPath')
    get_goog_trends(cells.value.split(','), cells_toCSV_FileName.value)

    cells_toCSV_FileName, ws = get_value_by_named_range(path, 'TX_FullPath')
    get_tx_fees(cells_toCSV_FileName.value)

    cells_toCSV_FileName, ws = get_value_by_named_range(path, 'BNB_Futures_FullPath')
    get_futures_binance(cells_toCSV_FileName.value)

    t1 = time.time()
    print('Entire crypto data collection took {} seconds'.format(t1-t0))  
