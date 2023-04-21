import pandas as pd
import option_config as cfg
from bs4 import BeautifulSoup
from selenium import webdriver
import openpyxl
import time

def get_value_by_named_range(strPath, strNamedRange):
    """Get Provided Named Range Cell From Given Excel File
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

class Ivolatility:
    Username = cfg.website['ivolatility']['username']
    Password = cfg.website['ivolatility']['password']
    strRootURL = 'https://www.ivolatility.com'
    driver = None

    def __init__(self, strSymbol, strExchange):
        self.symbol = strSymbol
        self.exchange = strExchange
        self.volatility = None

    @classmethod
    def signIn(cls):
        strURL = f'{cls.strRootURL}/login.j' 
        # Must Have Latest Stable Version Of Chrome Driver https://chromedriver.chromium.org/downloads
        cls.driver = webdriver.Chrome('Driver Location')
        cls.driver.get(strURL)
        time.sleep(10)
        username = cls.driver.find_element_by_name("username")
        username.send_keys(cls.Username)

        password = cls.driver.find_element_by_name("password")
        password.send_keys(cls.Password)

        cls.driver.find_element_by_class_name("btn-login").click()

    @classmethod
    def signOut(cls):
        cls.driver.quit()

    def getVolatility(self):
        strURL = f'{self.strRootURL}/options.j?ticker={self.symbol}:{self.exchange}&R=1'
        self.driver.get(strURL)
        time.sleep(10) # Wait 10 seconds for frame and ajax elements to load
        soup = BeautifulSoup(self.driver.page_source,'html.parser')
        tables = soup.findAll("table",{"class": "table-data"})
        data_list = pd.read_html(tables[1].prettify())
        dfVolatility = pd.concat(data_list, axis=1)
        dfVolatility.columns = dfVolatility.iloc[0]
        dfVolatility.insert(0, 'Ticker', self.symbol)
        self.volatility = dfVolatility


if __name__ == '__main__':
    # Collect Tickers From Portfolio Monitor (Excel)
    t0 = time.time()
    path = r'Path Location'
    
    cells_toCSV_FileName, ws = get_value_by_named_range(path, 'IVoption_Summary_FullPath')
    cells, ws = get_value_by_named_range(path, 'IV_TickerList')
    ticker_list = cells.value.split(",")

    cells, ws = get_value_by_named_range(path, 'IV_ExchangeList')
    exchange_list = cells.value.split(",")

    zipped_list = list(zip(ticker_list,exchange_list))

    Ivolatility.signIn()    

    optionData = []
    for ticker, exchange in zipped_list:
        option = Ivolatility(ticker,exchange)
        option.getVolatility()
        print('Successful import of {} data'.format(ticker))
        optionData.append(option.volatility)

    #Create Data Frame with data and save off as CSV
    df_optionData = pd.concat(optionData)
    df_optionData.to_csv(cells_toCSV_FileName.value)

    Ivolatility.signOut() 

    t1 = time.time()
    print('Entire IV option data collection took {} seconds'.format(t1-t0)) 