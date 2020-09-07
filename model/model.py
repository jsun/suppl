import os
import numpy as np
import pandas as pd
import torch
import torch.optim
import torch.nn
import torch.nn.functional as F
import torchvision
import torch.autograd

# os.environ['KMP_DUPLICATE_LIB_OK']='True'


PREFECTURE_DICT = {
                   'Hokkaido': 1,
                   'Aomori' : 2,
                   'Iwate' : 3,
                   'Miyagi' : 4,
                   'Akita' : 5,
                   'Yamagata' : 6,
                   'Fukushima': 7,
                   'Ibaraki': 8,
                   'Tochigi': 9,
                   'Gunma': 10,
                   'Saitama': 11,
                   'Chiba': 12,
                   'Tokyo': 13,
                   'Kanagawa': 14,
                   'Niigata': 15,
                   'Toyama': 16,
                   'Ishikawa': 17,
                   'Fukui': 18,
                   'Yamanashi': 19,
                   'Nagano': 20,
                   'Gifu': 21,
                   'Shizuoka': 22,
                   'Aichi': 23,
                   'Mie': 24,
                   'Shiga': 25,
                   'Kyoto': 26,
                   'Osaka': 27,
                   'Hyogo': 28,
                   'Nara': 29,
                   'Wakayama': 30,
                   'Tottori': 31,
                   'Shimane': 32,
                   'Okayama': 33,
                   'Hiroshima': 34,
                   'Yamaguchi': 35,
                   'Tokushima': 36,
                   'Kagawa': 37,
                   'Ehime': 38,
                   'Kochi': 39,
                   'Fukuoka': 40,
                   'Saga': 41,
                   'Nagasaki': 42,
                   'Kumamoto': 43, 
                   'Oita': 44,
                   'Miyazaki': 45,
                   'Kagoshima': 46,
                   'Okinawa': 47
                  }







class JppnetArch(torch.nn.Module):
    
    def __init__(self, n_hidden=12, dropout=0, activate_func='relu'):
        super(JppnetArch, self).__init__()
        
        if isinstance(n_hidden, int):
            n_hidden = [n_hidden]
        
        # input layer
        self.input = torch.nn.Linear(16, n_hidden[0])

        # hidden layer
        hidden_layers = []
        if len(n_hidden) > 1:
            for i in range(1, len(n_hidden)):
                if activate_func == 'relu':
                    hidden_layers += [torch.nn.ReLU(inplace=True)]
                else:
                    hidden_layers += [torch.nn.Sigmoid()]
                if dropout > 0:
                    hidden_layers += [torch.nn.Dropout(p=dropout)]
                hidden_layers += [torch.nn.Linear(n_hidden[i - 1], n_hidden[i])]

        if activate_func == 'relu':
            hidden_layers += [torch.nn.ReLU(inplace=True)]
        else:
            hidden_layers += [torch.nn.Sigmoid()]
        if dropout > 0:
            hidden_layers += [torch.nn.Dropout(p=dropout)]
        self.fc = torch.nn.Sequential(*hidden_layers)

        # output layer
        self.output = torch.nn.Linear(n_hidden[-1], 1)
    
    
    def forward(self, x):
        x = self.input(x)
        x = self.fc(x)
        x = self.output(x)
        return x






class textDataset(torch.utils.data.Dataset):
    
    
    def __init__(self, dat_fpath, is_test=False):
        self.is_test = is_test
        self.X, self.y = self.__read_csv(dat_fpath)
        self.size = len(self.y)
        
    
    def __read_csv(self, dat_fpath):
        '''
        load dataset from CSV file
        '''
        
        X = []
        y = []
        
        with open(dat_fpath, 'r') as datfh:
            for dat_buf in datfh:
                dat_buf = dat_buf.replace('\n', '').split('\t')
                if dat_buf[0] == 'incidence':
                    continue
                
                _y, year, month, prefecture, temp, rain, lat, lng = dat_buf
                
                # change month and prefecture name into one-hot codes
                _x = []
                _x.extend(self.__month2onehot(month))
                #_x.extend(self.__prefecture2onehot(prefecture))
                _x.extend([float(temp), float(rain), float(lat), float(lng)])
                
                if self.is_test:
                    X.append(_x)
                    y.append('NA')
                else:
                    if _y != 'NA':
                        _y = [float(_y)]
                        #_y = [np.log10(float(_y) + 1)]
                        X.append(_x)
                        y.append(_y)
        
        X = torch.from_numpy(np.array(X)).float()
        y = torch.from_numpy(np.array(y)).float()
        
        return X, y
        
    
    def __len__(self):
        return self.y.shape[0]
    
    
    def __getitem__(self, i):
        _x = self.X[i]
        _y = self.y[i]
        
        if self.is_test:
            return _x
        else:
            return _x, _y
    
    
    def __month2onehot(self, m):
        m_vec = np.zeros(12)
        m_vec[int(m) - 1] = 1
        return m_vec.tolist()
        #m = int(m)
        #if m == 12 or m == 1 or m == 2:
        #    m_vec = [1, 0, 0, 0]
        #elif m == 3 or m == 4 or m == 5:
        #    m_vec = [0, 1, 0, 0]
        #elif m == 6 or m == 7 or m == 8:
        #    m_vec = [0, 0, 1, 0]
        #elif m == 9 or m == 10 or m == 11:
        #    m_vec = [0, 0, 0, 1]
        #elif m == 9 or m == 10 or m == 11:
        #    m_vec = [0, 0, 0, 1]
        #elif m == 9 or m == 10 or m == 11:
        #    m_vec = [0, 0, 0, 1]
        #
        #return m_vec
    
    def __prefecture2onehot(self, p):
        p_vec = np.zeros(len(PREFECTURE_DICT))
        p_vec[PREFECTURE_DICT[p] - 1] = 1
        return p_vec.tolist()







class Jppnet():
    
    def __init__(self, n_hidden=12, dropout=0.5, activate_func='relu'):
        
        self.model = JppnetArch(n_hidden=n_hidden, dropout=dropout, activate_func=activate_func) 
        self.__best_model = JppnetArch(n_hidden=n_hidden, dropout=dropout, activate_func=activate_func) 
        #self.device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
        self.device = torch.device('cpu')
        self.model.to(self.device)
        
        
    def load_weights(self, weights_fpath):
        self.model.load_state_dict(torch.load(weights_fpath))
    
    
    def save(self, weights_fpath):
        torch.save(self.model.state_dict(), weights_fpath)
    
    
    def train(self, train_data, valid_data, batch_size=1024, num_epochs=100, num_workers=4):
        
        # load train and validation dataset
        datasets = {
            'train': textDataset(train_data),
            'valid': textDataset(valid_data)
        }
        dataloader = {
            'train': torch.utils.data.DataLoader(datasets['train'],
                                                  batch_size=batch_size, shuffle=True, num_workers=num_workers),
            'valid': torch.utils.data.DataLoader(datasets['valid'],
                                                  batch_size=batch_size, num_workers=num_workers)
        }
        dataset_sizes = {
            'train': len(datasets['train']),
            'valid': len(datasets['valid']),
        }
        
        
        # training parameters
        criterion = torch.nn.MSELoss()
        optimizer = torch.optim.Adam(self.model.parameters(), lr=1e-2)
        scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=20, gamma=0.1)
        
        # training and validation
        epoch_loss = {'train': [], 'valid': []}
        _last_epoch_loss = 1e10
        _last_epoch_loss_count = 0
        for epoch in range(num_epochs):
            
            for phase in ['train', 'valid']:

                if phase == 'train':
                    self.model.train()
                else:
                    self.model.eval()
            
                running_loss = 0.0
                running_corrects = 0
        
                for inputs, labels in dataloader[phase]:
                    inputs = inputs.to(self.device)
                    labels = labels.to(self.device)
                    
                    # forward
                    with torch.set_grad_enabled(phase == 'train'):
                        outputs = self.model(inputs)
                        loss = criterion(outputs, labels)
                        
                        if phase == 'train':
                            optimizer.zero_grad()
                            loss.backward()
                            optimizer.step()
                
                    running_loss += loss.item() * inputs.size(0)
                
                if phase == 'train':
                    scheduler.step()

                _epoch_loss = running_loss / dataset_sizes[phase]
                
                if phase == 'valid':
                    if _epoch_loss < _last_epoch_loss:
                        _last_epoch_loss = _epoch_loss
                        _last_epoch_loss_count = 0
                        self.__best_model.load_state_dict(self.model.state_dict())
                    else:
                        _last_epoch_loss_count += 1
            
                #print('{} Loss: {:.4f}'.format(phase, _epoch_loss))
                epoch_loss[phase].append(_epoch_loss)
            
            if _last_epoch_loss_count > 20:
                break
        
        epoch_loss = pd.DataFrame(epoch_loss)
        self.model.load_state_dict(self.__best_model.state_dict())
        return epoch_loss
    
    
    
    
    def validate(self, test_data, batch_size=None):
        
        self.model.eval()
        
        dataset = textDataset(test_data)
        dataloader = torch.utils.data.DataLoader(dataset, batch_size=1024, shuffle=False)
        
        pred = {'label': [], 'predicted': []}
        
        with torch.set_grad_enabled(False):
            for inputs, labels in dataloader:
                inputs = inputs.to(self.device)
                pred['predicted'].extend(self.model(inputs).data.cpu().numpy().flatten())
                pred['label'].extend(labels.data.cpu().numpy().flatten())
        
        pred = pd.DataFrame(pred)
        
        return pred
    
    
    
     def inference(self, test_data, batch_size=None):
        
        self.model.eval()
        
        dataset = textDataset(test_data, is_test)
        dataloader = torch.utils.data.DataLoader(dataset, batch_size=1024, shuffle=False)
        
        pred = {'label': [], 'predicted': []}
        
        with torch.set_grad_enabled(False):
            for inputs, labels in dataloader:
                inputs = inputs.to(self.device)
                pred['predicted'].extend(self.model(inputs).data.cpu().numpy().flatten())
                pred['label'].extend(labels.data.cpu().numpy().flatten())
        
        pred = pd.DataFrame(pred)
        
        return pred
    
    
    
    

