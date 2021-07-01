import numpy as np
import os
from shutil import copy2
from sklearn import preprocessing
import tensorflow as tf

def split_dataset(BASE_PATH = 'flower_photos', DATASET_PATH = 'dataset', train_images = 0.8, val_images = 0.1, test_images = 0.1):
    # Specify path to the downloaded folder
    classes = os.listdir(BASE_PATH)

    # Specify path for copying the dataset into train and val sets
    if os.path.exists(DATASET_PATH):
        os.system("rm -rf "+DATASET_PATH)
    os.makedirs(DATASET_PATH, exist_ok=True)

    # Creating train directory
    train_dir = os.path.join(DATASET_PATH, 'train')
    os.makedirs(train_dir, exist_ok=True)

    # Creating val directory
    val_dir = os.path.join(DATASET_PATH, 'val')
    os.makedirs(val_dir, exist_ok=True)

    # Creating test directory
    test_dir = os.path.join(DATASET_PATH, 'test')
    os.makedirs(test_dir, exist_ok=True)   

    # Copying images from original folder to dataset folder
    for class_name in classes:
        if len(class_name.split('.')) >= 2:
            continue
        print(f"Copying images for {class_name}...")
        

        # Creating destination folder (train and val)
        class_train_dir = os.path.join(train_dir, class_name)
        os.makedirs(class_train_dir, exist_ok=True)
        
        class_val_dir = os.path.join(val_dir, class_name)
        os.makedirs(class_val_dir, exist_ok=True)

        class_test_dir = os.path.join(test_dir, class_name)
        os.makedirs(class_test_dir, exist_ok=True)

        # Shuffling the image list
        class_path = os.path.join(BASE_PATH, class_name)
        class_images = os.listdir(class_path)
        tot_file = 0
        for Files in class_images:
            tot_file += 1

        np.random.shuffle(class_images)

        for image in class_images[:int(tot_file*train_images)]:
            copy2(os.path.join(class_path, image), class_train_dir)
        for image in class_images[int(tot_file*train_images):int(tot_file*train_images)+int(tot_file*val_images)]:
            copy2(os.path.join(class_path, image), class_val_dir)
        for image in class_images[int(tot_file*train_images)+int(tot_file*val_images):]:
            copy2(os.path.join(class_path, image), class_test_dir)


class Generator(tf.keras.utils.Sequence):

    def __init__(self, DATASET_PATH, BATCH_SIZE=32, shuffle_images=True, image_min_side=8):
        """ Initialize Generator object.
        Args
            DATASET_PATH           : Path to folder containing individual folders named by their class names
            BATCH_SIZE             : The size of the batches to generate.
            shuffle_images         : If True, shuffles the images read from the DATASET_PATH
            image_min_side         : After resizing the minimum side of an image is equal to image_min_side.
        """

        self.batch_size = BATCH_SIZE
        self.shuffle_images = shuffle_images
        self.image_min_side = image_min_side
        self.load_image_paths_labels(DATASET_PATH)
        self.create_image_groups()
    
    def load_image_paths_labels(self, DATASET_PATH):
        
        classes = os.listdir(DATASET_PATH)
        lb = preprocessing.LabelBinarizer()
        lb.fit(classes)

        self.image_paths = []
        self.image_labels = []
        for class_name in classes:
            class_path = os.path.join(DATASET_PATH, class_name)
            for image_file_name in os.listdir(class_path):
                self.image_paths.append(os.path.join(class_path, image_file_name))
                self.image_labels.append(class_name)

        self.image_labels = np.array(lb.transform(self.image_labels), dtype='float32')
        
        assert len(self.image_paths) == len(self.image_labels)

    def create_image_groups(self):
        if self.shuffle_images:
            #Randomly shuffle dataset
            seed = 4321
            np.random.seed(seed)
            np.random.shuffle(self.image_paths)
            np.random.seed(seed)
            np.random.shuffle(self.image_labels)

        # Divide image_paths and image_labels into groups of BATCH_SIZE
        self.image_groups = [[self.image_paths[x % len(self.image_paths)] for x in range(i, i + self.batch_size)]
                              for i in range(0, len(self.image_paths), self.batch_size)]
        self.label_groups = [[self.image_labels[x % len(self.image_labels)] for x in range(i, i + self.batch_size)]
                              for i in range(0, len(self.image_labels), self.batch_size)]


    def load_images(self, image_group):

        images = []
        for image_path in image_group:
            img = np.genfromtxt(image_path,delimiter=',',skip_header=1)
            img = np.absolute(img)
            img = img[..., np.newaxis]
            images.append(img)

        return images

    def construct_image_batch(self, image_group):
        #max_shape = tuple(max(image.shape[x] for image in image_group) for x in range(3))
        max_shape = (9,200,1)
        # construct an image batch object
        image_batch = np.zeros((self.batch_size,) + max_shape, dtype='float32')

        # copy all images to the upper left part of the image batch object
        for image_index, image in enumerate(image_group):
            image_batch[image_index, :image.shape[0], :image.shape[1], :image.shape[2]] = image

        return image_batch
    
    def __len__(self):
        """
        Number of batches for generator.
        """

        return len(self.image_groups)

    def __getitem__(self, index):
        """
        Keras sequence method for generating batches.
        """
        image_group = self.image_groups[index]
        label_group = self.label_groups[index]
        images = self.load_images(image_group)
        image_batch = self.construct_image_batch(images)

        return np.array(image_batch), np.array(label_group)


model= tf.keras.models.Sequential()
model.add(tf.keras.layers.Input(shape=(9, 200, 1)))


model.add(tf.keras.layers.Flatten())


model.add(tf.keras.layers.Dense(32))
model.add(tf.keras.layers.Dropout(rate=0.2))
model.add(tf.keras.layers.BatchNormalization())
model.add(tf.keras.layers.Activation('relu'))


model.add(tf.keras.layers.Dense(1))

model.add(tf.keras.layers.Activation('sigmoid'))


callback = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=3,restore_best_weights=True)


model.compile(loss='binary_crossentropy',
optimizer=tf.keras.optimizers.SGD(lr=0.01,decay=1e-6, momentum=0.9,nesterov=True),
metrics='binary_accuracy')

split_dataset('Train_dataset',DATASET_PATH='dataset',train_images = 0.8, val_images = 0.1,test_images = 0.1)


train_dir = 'dataset/train'
val_dir = 'dataset/val'

BATCH_SIZE=80


train_generator = Generator(train_dir, BATCH_SIZE, shuffle_images=True, image_min_side=1)
val_generator = Generator(val_dir, 80, shuffle_images=True, image_min_side=1)


history = model.fit(train_generator,
steps_per_epoch=len(train_generator),
epochs=50,
verbose=1,
callbacks=[callback],
validation_data = val_generator,
validation_steps=len(val_generator))

model.save('model_smRNA')


