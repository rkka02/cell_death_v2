import torch
from tqdm import tqdm
import datetime

def train_model(model, train_loader, val_loader, criterion, optimizer, 
                num_epochs:int=20, save_model:bool=False):
    # Train
    num_epochs = num_epochs

    for epoch in tqdm(range(num_epochs)):
        # train
        model.train()
        train_loss, train_correct, train_total = 0, 0, 0
        # i=1
        for images, labels in train_loader:
            # print(f'batch {i} started! {datetime.datetime.now()}')
            # i += 1
            images, labels = images.cuda(), labels.cuda()
            outputs = model(images)
            loss = criterion(outputs, labels)
            train_loss += loss.item()
            
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            
            _, preds = torch.max(outputs, 1)
            train_correct += (preds==labels).sum().item()
            train_total += len(labels)
        
        #validation
        model.eval()
        val_loss, val_correct, val_total = 0, 0, 0
        with torch.no_grad():
            for images, labels in val_loader:
                images, labels = images.cuda(), labels.cuda()
                
                outputs = model(images)
                loss = criterion(outputs, labels)
                val_loss += loss
                
                _, preds = torch.max(outputs, 1)
                val_correct += (preds==labels).sum().item()
                val_total += len(labels)

        # Save model if save_model==True
        if save_model==True:        
            torch.save(model.state_dict(), f'epoch_{epoch}_3d_val_acc_{val_correct/val_total:.4f}.pth') 
              
        print(f"Epoch : {epoch}")
        print(f"train loss : {train_loss/train_total:.6f} || train_acc : {train_correct/train_total:.4f}")
        print(f"val loss : {val_loss/val_total:.6f} || val_acc : {val_correct/val_total:.4f}")
        