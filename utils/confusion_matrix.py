import torch
import seaborn as sns
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
import matplotlib.pyplot as plt

def confusion_matrix(model, val_loader):
    model = model.cuda()

    model.eval()
    all_preds = []
    all_labels = []

    with torch.no_grad():
        for images, labels in val_loader:
            images, labels = images.cuda(), labels.cuda()

            outputs = model(images)
            _, preds = torch.max(outputs, 1)

            all_preds.extend(preds.cpu().numpy())
            all_labels.extend(labels.cpu().numpy())

    # Calculate confusion matrix
    cm = confusion_matrix(all_labels, all_preds, normalize='true')

    # Plot confusion matrix
    plt.figure(figsize=(8, 6))
    sns.heatmap(cm, annot=True, cmap="Blues", xticklabels=["Apoptosis", "Necroptosis", "Necrosis", "Live"], 
                                            yticklabels=["Apoptosis", "Necroptosis", "Necrosis", "Live"], annot_kws={'size': 24})
    plt.xlabel("Predicted", fontsize=24, labelpad=20)
    plt.ylabel("Ground Truth", fontsize=24, labelpad=20)
    plt.title("Confusion Matrix", fontsize=24, pad=20)

    # Adjust colorbar font size
    cbar = plt.gca().collections[0].colorbar
    cbar.ax.tick_params(labelsize=18)  # Set colorbar tick font size

    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    plt.savefig('figures/figure3/confusion_matrix.svg', format='svg', transparent=True)
    plt.show()

    print('Confusion matrix is successfully generated')